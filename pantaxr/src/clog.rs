
use flexi_logger::{DeferredNow, Record, style};

// chrono only support %.3f
const TS_DASHES_BLANK_COLONS_DOT_BLANK: &str = "%Y-%m-%d %H:%M:%S%.3f";

pub fn init_logger(level: &str) {
    flexi_logger::Logger::try_with_str(level)
        .expect("Logger config failed")
        .log_to_stdout()
        .set_palette("196;208;51;129;8".to_string())
        .format(my_own_format_colored)
        .start()
        .expect("Logger start failed");
}

fn my_own_format_colored(
    w: &mut dyn std::io::Write,
    now: &mut DeferredNow,
    record: &Record,
) -> Result<(), std::io::Error> {
    let paintlevel = record.level();

    write!(
        w,
        "{} {} [{}] {}",
        now.format(TS_DASHES_BLANK_COLONS_DOT_BLANK),
        style(paintlevel).paint(record.level().to_string()),
        record.module_path().unwrap_or(""),
        &record.args()
    )
}