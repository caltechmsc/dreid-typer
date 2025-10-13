mod harness;

use harness::cases::amino_acids::*;
use harness::run_molecule_test_case;

macro_rules! generate_molecule_test {
    ($test_name:ident, $molecule_case:expr) => {
        #[test]
        fn $test_name() {
            run_molecule_test_case(&$molecule_case);
        }
    };
}

generate_molecule_test!(glycine_zwitterion_is_typed_correctly, GLYCINE_ZWITTERION);
generate_molecule_test!(alanine_zwitterion_is_typed_correctly, ALANINE_ZWITTERION);
