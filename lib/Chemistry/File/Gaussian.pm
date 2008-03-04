package Chemistry::File::Gaussian;

$VERSION = '0.10';
# $Id$

use 5.006;
use strict;
use warnings;
use base "Chemistry::File";
use Chemistry::Mol 0.25;
use Chemistry::InternalCoords;
use Carp;

=head1 NAME

Chemistry::File::Gaussian - Gaussian output file reader

=head1 SYNOPSIS

    use Chemistry::File::Gaussian;

    # read an output file
    my $mol = Chemistry::Mol->read('file.out');

    # write an input file using cartesian coordinates
    $mol->write('file.com');

=cut

=head1 DESCRIPTION

=cut

Chemistry::Mol->register_format("g03");

sub slurp_mol {
    my ($self, $fh, %opts) = @_;
    return if $fh->eof;
    my $s;
    if ($self->{gaussian_type} eq 'input') {
        $s = $self->slurp;
    } else {
        if (@{$self->mols} == 0) {
            while (<$fh>) {
                last if /^$/;
                $s .= $_;
            }
        } else {
            while (<$fh>) {
                last if /orientation/;
            }
            <$fh> for 1 .. 4;  # skip table header
            while (<$fh>) {
                last if /-----------/;
                $s .= $_;
            }
        }
    }
    return $s;
}

sub read_mol {
    my ($self, $fh, %opts) = @_;

    return if $fh->eof;
    my $mol_class  = $opts{mol_class}  || "Chemistry::Mol";
    my $atom_class = $opts{atom_class} || $mol_class->atom_class;
    my $bond_class = $opts{bond_class} || $mol_class->bond_class;
    my $mol;

    if ($self->{gaussian_type} eq 'input' or @{$self->mols} == 0) {
        $mol = $mol_class->new();
        my %vars;
        my @var_list;
        # read zmat...
        while (<$fh>) {
            last if /(?:Charge =)? *([0-9+-]+) +(?:Multiplicity =)? *(\d+)/;
        }
        my ($charge, $multiplicity) = ($1, $2);

        # read zmat block
        while (<$fh>) {
            last if /^ *$/ or /Variables:/;
            my ($sym, $len_ref, $len, $ang_ref, $ang, $dih_ref, $dih) = split;
            $sym = ucfirst(lc($sym));
            my $atom = $mol->new_atom(symbol => $sym);
            for my $var ($len, $ang, $dih) {
                next unless defined $var;
                if ($var =~ /[A-Za-z]/) {
                    push @var_list, $var unless exists $vars{$var};
                    $vars{$var} = undef;
                }
            }
            $atom->internal_coords($len_ref, $len, $ang_ref, $ang, 
                $dih_ref, $dih);
        }

        $mol->attr('gaussian/var_list', \@var_list);
        $mol->attr('gaussian/vars', \%vars);

        # read variables
        while (<$fh>) {
            last if /^ *$/;
            next if /Constants:/; # XXX Treat constants like variables
            my ($name, $val) = split;
            $vars{$name} = $val;
        }

        # set variables and generate cartesians
        for my $atom ($mol->atoms) {
            my $ic = $atom->internal_coords;
            for my $type (qw(distance angle dihedral)) {
                my ($ref, $name) = $ic->$type;
                next unless defined $name and exists $vars{$name};
                $ic->$type($vars{$name});
                $ic->attr("$type name", $name);
            }
            $ic->add_cartesians;
        }

        if ($self->{gaussian_type} eq 'input') {
            local $/; <$fh> # slurp
        }
    } else { # non-first molecules in output file
        $mol = $self->mols->[0]->clone;

        while (<$fh>) {
            last if /orientation/;
        }
        return if $fh->eof;
        <$fh> for 1 .. 4;  # skip table header
        my $i = 0;
        while (<$fh>) {
            last if /-----------/;
            $i++;
            my (undef,undef,undef,@coords) = split;
            my $atom = $mol->atoms($i);
            $atom->coords(@coords);
            my $ic = $atom->internal_coords;
            $ic->update;

            for my $type (qw(distance angle dihedral)) {
                my ($ref, $val) = $ic->$type;
                my $name = $ic->attr("$type name");
                my $vars = $mol->attr("gaussian/vars");
                next unless defined $name and exists $vars->{$name};
                $vars->{$name} = $val;
            }

        }
    }
    return $mol;
}


sub read_header {
    my ($self, %opts) = @_;
    my $fh = $self->fh;
    my $line = <$fh>;
    if ($line =~ /Entering Gaussian System/) {
        $self->{gaussian_type} = 'output';
        while (<$fh>) {
            last if /Symbolic Z-matrix/;
        }
    } else {
        $self->{gaussian_type} = 'input';
        seek $fh, 0, 0;
    }
}


sub file_is {
    my ($class, $fname) = @_;
    
    return 1 if $fname =~ /\.(?:com)$/i;

    open F, $fname or croak "Could not open file $fname";
    
    my $line = <F>;
    close F;
    return 1 if $line =~ /Entering Gaussian System/;
    return 0;
}

sub name_is {
    my ($class, $fname) = @_;
    $fname =~ /\.(?:out|com)$/i;
}

sub write_string {
    my ($class, $mol, %opts) = @_;
    %opts = (coords => "internal", rebuild => 0, %opts);
    my $ret = $class->format_header($mol);
    if ($opts{coords} eq 'cartesian') {
        $ret .= $class->format_line_cart($_) for $mol->atoms;
    } else {
        if ($opts{rebuild} or ! $mol->atoms(1)->internal_coords ) {
            require Chemistry::InternalCoords::Builder;
            Chemistry::InternalCoords::Builder::build_zmat($mol);
        }
        my %index;
        @index{$mol->atoms} = (1 .. $mol->atoms);
        $class->name_vars($mol) if $opts{name_vars};
        $ret .= $class->format_line_ic($_, \%index) for $mol->atoms;
        $ret .= $class->format_vars($mol);
    }
    $ret;
}


sub format_header {
    my ($class, $mol) = @_;
    my $ret = ($mol->attr("gaussian/keywords") || '#') . "\n\n";
    my $name = $mol->name || '';
    if ($name) {
        $name =~ s/\n$/ /s;
        $ret .= $name . "\n\n";
    }
    $ret .= $mol->charge . " " . ($mol->attr('multiplicity') || 1) . "\n";
    $ret;
}

sub format_line_ic {
    my ($class, $atom, $index) = @_;
    my $ic = $atom->internal_coords;
    my ($len_ref, $len_val, $len_name) = ($ic->distance, 
        $ic->attr('distance name'));
    my ($ang_ref, $ang_val, $ang_name) = ($ic->angle,    
        $ic->attr('angle name'));
    my ($dih_ref, $dih_val, $dih_name) = ($ic->dihedral, 
        $ic->attr('dihedral name'));
    my $len_idx = $index->{$len_ref||0} || 0;
    my $ang_idx = $index->{$ang_ref||0} || 0;
    my $dih_idx = $index->{$dih_ref||0} || 0;
    no warnings 'uninitialized';

    my $line = sprintf "%-2s", $atom->symbol;
    
    if ($len_idx) {
        $line .= sprintf " %4d", $len_idx;
        $line .= $len_name ? sprintf " %-8s", $len_name 
            : sprintf " %8.3f", $len_val;
    }

    if ($ang_idx) {
        $line .= sprintf " %4d", $ang_idx;
        $line .= $ang_name ? sprintf " %-8s", $ang_name 
            : sprintf " %8.3f", $ang_val;
    }
    
    if ($dih_idx) {
        $line .= sprintf " %4d", $dih_idx;
        $line .= $dih_name ? sprintf " %-8s", $dih_name 
            : sprintf " %8.3f", $dih_val;
    }

    $line .= "\n";

    return $line;
}

sub format_line_cart {
    my ($class, $atom) = @_;
    my ($x, $y, $z) = $atom->coords->array;
    sprintf "%-2s %8.3f %8.3f %8.3f \n",
        $atom->symbol, 
        $atom->coords->array;
}

sub name_vars {
    my ($self, $mol) = @_;

    my %vars;
    my @var_list;
    my $i;
    my %prefix = (distance => 'b', angle => 'a', dihedral => 'd');
    for my $atom ($mol->atoms) {
        $i++;
        my $ic = $atom->internal_coords;
        for my $type (qw(distance angle dihedral)) {
            my ($ref, $val) = $ic->$type; 
            next unless $ref;
            my $name = $prefix{$type} . $i;
            $ic->attr("$type name", $name);
            push @var_list, $name;
            $vars{$name} = $val;
        }
    }
    $mol->attr('gaussian/vars', \%vars);
    $mol->attr('gaussian/var_list', \@var_list);
}

sub format_vars {
    my ($class, $mol) = @_;
    my $vars   = $mol->attr('gaussian/vars');
    my $ret = "\n";
    if ($vars) {
        for my $name (@{$mol->attr('gaussian/var_list')}){
            $ret .= sprintf "%-8s %8.3f\n", $name, $vars->{$name};
        }
    }
    $ret .= "\n\n";
    # XXX TODO : constants
    #my $consts = $mol->attr('gaussian/consts');
    return $ret;
}

1;

=head1 TO DO

When writing a Mopac file, this version marks all coordinates as variable 
(for the purpose of geometry optimization by Mopac). A future version should
have more flexibility.

=head1 VERSION

0.15

=head1 SEE ALSO

L<Chemistry::Mol>, L<Chemistry::File>, L<Chemistry::InternalCoords>, 
L<Chemistry::InternalCoords::Builder>, L<http://www.perlmol.org/>.

=head1 AUTHOR

Ivan Tubert-Brohman <itub@cpan.org>

=head1 COPYRIGHT

Copyright (c) 2004 Ivan Tubert. All rights reserved. This program is free
software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut

