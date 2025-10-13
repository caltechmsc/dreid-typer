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
generate_molecule_test!(valine_zwitterion_is_typed_correctly, VALINE_ZWITTERION);
generate_molecule_test!(leucine_zwitterion_is_typed_correctly, LEUCINE_ZWITTERION);
generate_molecule_test!(
    isoleucine_zwitterion_is_typed_correctly,
    ISOLEUCINE_ZWITTERION
);
generate_molecule_test!(proline_zwitterion_is_typed_correctly, PROLINE_ZWITTERION);
generate_molecule_test!(serine_zwitterion_is_typed_correctly, SERINE_ZWITTERION);
generate_molecule_test!(
    threonine_zwitterion_is_typed_correctly,
    THREONINE_ZWITTERION
);
generate_molecule_test!(cysteine_zwitterion_is_typed_correctly, CYSTEINE_ZWITTERION);
generate_molecule_test!(
    methionine_zwitterion_is_typed_correctly,
    METHIONINE_ZWITTERION
);
generate_molecule_test!(
    aspartate_zwitterion_is_typed_correctly,
    ASPARTATE_ZWITTERION
);
generate_molecule_test!(
    asparagine_zwitterion_is_typed_correctly,
    ASPARAGINE_ZWITTERION
);
generate_molecule_test!(
    glutamate_zwitterion_is_typed_correctly,
    GLUTAMATE_ZWITTERION
);
generate_molecule_test!(
    glutamine_zwitterion_is_typed_correctly,
    GLUTAMINE_ZWITTERION
);
generate_molecule_test!(lysine_zwitterion_is_typed_correctly, LYSINE_ZWITTERION);
generate_molecule_test!(arginine_zwitterion_is_typed_correctly, ARGININE_ZWITTERION);
