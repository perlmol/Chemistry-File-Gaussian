#!/usr/bin/perl

use warnings;
use strict;
use Chemistry::File qw(MDLMol XYZ PDB);
use Chemistry::Bond::Find qw(find_bonds);
use lib 'lib';
use Chemistry::File::Gaussian;
use diagnostics;

#my $mol = Chemistry::Mol->read('cyclo-sto.com');

#print $mol->print(format => "g03");
#print $mol->print(format => "xyz");

local $SIG{__WARN__} = sub { die $_[0] };

my @mols = Chemistry::Mol->read(shift || 'cyclo-sto.out');

print "got ". @mols . " structures\n";

my $i = 0;
for my $mol (@mols) {
    $i++;
    my $name = sprintf "g%02d", $i;
    $mol->write("$name.com", format => "g03");
    $mol->write("$name.pdb", format => "pdb");
}
