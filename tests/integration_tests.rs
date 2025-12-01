mod harness;

use harness::cases::amino_acids::*;
use harness::cases::dreiding_paper::*;
use harness::cases::nucleic_acids::*;
use harness::run_molecule_test_case;

macro_rules! generate_molecule_test {
    ($test_name:ident, $molecule_case:expr) => {
        #[test]
        fn $test_name() {
            run_molecule_test_case(&$molecule_case);
        }
    };
}
