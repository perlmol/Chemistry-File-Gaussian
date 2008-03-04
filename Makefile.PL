use 5.006;
use ExtUtils::MakeMaker;

WriteMakefile(
    'NAME'              => 'Chemistry::File::Gaussian',
    'VERSION_FROM'      => 'lib/Chemistry/File/Gaussian.pm',
    'PREREQ_PM'         => {
        "Chemistry::Mol" => "0.25",
        "Chemistry::InternalCoords" => '0.15',
        "List::Util" => 0,
    },
);
