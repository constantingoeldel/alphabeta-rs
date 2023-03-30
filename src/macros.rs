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
macro_rules! assert_close {
    ($x:expr, $y:expr ) => {
        if ($x - $y).abs() > 1e-4 {
            panic!(
                "assertion failed: `abs(left - right) < {}`, (left: `{}`, right: `{}`)",
                1e-6, $x, $y
            );
        }
    };
}

#[macro_export]
macro_rules! assert_close_10_percent {
    ($x:expr, $y:expr ) => {
        if ($x - $y).abs() > 0.1 * $y {
            panic!(
                "assertion failed: `abs(left - right) < 10 % `, (left: `{}`, right: `{}`)",
                $x, $y
            );
        }
    };
}
