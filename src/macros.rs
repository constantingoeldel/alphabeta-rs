#[macro_export]
macro_rules! assert_close_epsilon {
    ($x:expr, $y:expr, $d:expr) => {
        if !(($x - $y).abs() < $d) {
            panic!(
                "assertion failed: `abs(left - right) < {}`, (left: `{}`, right: `{}`)",
                $d, $x, $y
            );
        }
    };
}

#[macro_export]
macro_rules! assert_within_10_percent {
    ($x:expr, $y:expr ) => {
        if ($x - $y).abs() > 0.1 * $y {
            panic!(
                "assertion failed: `abs(left - right) < 10 % `, (left: `{}`, right: `{}`)",
                $x, $y
            );
        }
    };
}
