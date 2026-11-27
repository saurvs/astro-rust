//! Meeus, chapter 9: Jewish and Moslem (Islamic) calendars.

use astro::calendar::{self, IslamicDate, JewishYearKind};
use astro::time::{self, CalType, Date, Month, Weekday};

fn greg(year: i16, month: Month, day: f64) -> Date {
    Date {
        year,
        month,
        decimal_day: day,
        cal_type: CalType::Gregorian,
    }
}

fn jul(year: i16, month: Month, day: f64) -> Date {
    Date {
        year,
        month,
        decimal_day: day,
        cal_type: CalType::Julian,
    }
}

/// Fixed-day (R.D.) number to JD at 0h, as in Reingold & Dershowitz,
/// *Calendrical Calculations*.
fn jd_from_rd(rd: i64) -> f64 {
    rd as f64 + 1721424.5
}

/// Reference Hebrew calendar: R.D. of 1 Tishri (Reingold & Dershowitz).
fn hebrew_new_year_rd(year: i64) -> i64 {
    fn elapsed(y: i64) -> i64 {
        let months = (235 * y - 234).div_euclid(19);
        let parts = 12084 + 13753 * months;
        let day = 29 * months + parts.div_euclid(25920);
        if (3 * (day + 1)).rem_euclid(7) < 3 {
            day + 1
        } else {
            day
        }
    }
    let (ny0, ny1, ny2) = (elapsed(year - 1), elapsed(year), elapsed(year + 1));
    let correction = if ny2 - ny1 == 356 {
        2
    } else if ny1 - ny0 == 382 {
        1
    } else {
        0
    };
    -1373427 + elapsed(year) + correction
}

/// Reference arithmetical Islamic calendar (Reingold & Dershowitz).
const ISLAMIC_EPOCH_RD: i64 = 227015;

fn rd_from_islamic(y: i64, m: i64, d: i64) -> i64 {
    d + 29 * (m - 1)
        + (6 * m - 1).div_euclid(11)
        + (y - 1) * 354
        + (3 + 11 * y).div_euclid(30)
        + ISLAMIC_EPOCH_RD
        - 1
}

fn islamic_from_rd(rd: i64) -> (i64, i64, i64) {
    let y = (30 * (rd - ISLAMIC_EPOCH_RD) + 10646).div_euclid(10631);
    let prior = rd - rd_from_islamic(y, 1, 1);
    let m = (11 * prior + 330).div_euclid(325);
    let d = rd - rd_from_islamic(y, m, 1) + 1;
    (y, m, d)
}

// Meeus, example 9.a
#[test]
fn example_9a_pesach_1990() {
    let p = calendar::pesach(1990);
    assert_eq!(p, greg(1990, Month::Apr, 10.0));

    // 1 Tishri 5751 is 163 days later, 1990 September 20
    let ny = calendar::jewish_new_year(5751).unwrap();
    assert_eq!(ny, greg(1990, Month::Sept, 20.0));
    assert_eq!(time::days_between(&p, &ny), 163.0);
    // 5751 is a common year; 5752 begins on 1991 September 9
    assert!(!calendar::is_jewish_leap_year(5751));
    assert_eq!(
        calendar::jewish_new_year(5752).unwrap(),
        greg(1991, Month::Sept, 9.0)
    );
    let y = calendar::jewish_year(5751).unwrap();
    assert_eq!(y.new_year, ny);
    assert_eq!(y.days, 354);
    assert!(!y.embolismic);
    assert_eq!(y.kind, JewishYearKind::Regular);
}

#[test]
fn pesach_calendar_switch() {
    assert_eq!(calendar::pesach(1582).cal_type, CalType::Julian);
    assert_eq!(calendar::pesach(1583).cal_type, CalType::Gregorian);
    // a year spanning the change of calendar has a normal length
    let y = calendar::jewish_year(1582 + 3761).unwrap();
    assert!([353, 354, 355, 383, 384, 385].contains(&y.days));
}

// Pesach agrees with the Hebrew calendar of Reingold & Dershowitz (through
// 1 Tishri of the next year, which is always 163 days later).
#[test]
fn pesach_matches_reference_hebrew_calendar() {
    for year in -3700_i16..=29000 {
        let jd = time::julian_day(&calendar::pesach(year));
        let reference = jd_from_rd(hebrew_new_year_rd(year as i64 + 3761) - 163);
        assert_eq!(jd, reference, "{year}");
    }
}

// "Lo BaDU Pesach": 15 Nisan never falls on Monday, Wednesday or Friday, and
// "Lo ADU Rosh": 1 Tishri never falls on Sunday, Wednesday or Friday.
#[test]
fn pesach_and_new_year_weekdays() {
    for year in 1..=9999_i16 {
        let p = calendar::pesach(year);
        let w = time::weekday_frm_date(&p);
        assert!(
            !matches!(w, Weekday::Monday | Weekday::Wednesday | Weekday::Friday),
            "{year}: {w:?}"
        );
        let ny = calendar::jewish_new_year(year + 3761).unwrap();
        let w = time::weekday_frm_date(&ny);
        assert!(
            !matches!(w, Weekday::Sunday | Weekday::Wednesday | Weekday::Friday),
            "{year}: {w:?}"
        );
        // In Gregorian years of the present era the new year falls in
        // September or October. (Not always in Julian years: 1 Tishri 5261
        // was 1500 August 25, Julian. The Jewish year is also about 0.0047
        // day longer than the Gregorian year, so the date drifts later.)
        if (1583..=2500).contains(&year) {
            assert!(ny.month == Month::Sept || ny.month == Month::Oct, "{ny:?}");
        }
    }
}

#[test]
fn jewish_year_lengths_and_cycle() {
    let mut total_days = 0_u64;
    let mut leap = 0;
    let first = 5701;
    let cycles = 50;
    for a in first..first + 19 * cycles {
        let y = calendar::jewish_year(a).unwrap();
        assert_eq!(y.embolismic, y.days > 380, "{a}: {y:?}");
        let expected_kind = match y.days % 10 {
            3 => JewishYearKind::Deficient,
            4 => JewishYearKind::Regular,
            _ => JewishYearKind::Complete,
        };
        assert_eq!(y.kind, expected_kind);
        assert_eq!(y.new_year, calendar::jewish_new_year(a).unwrap(), "{a}");
        total_days += y.days as u64;
        if y.embolismic {
            leap += 1;
        }
    }
    // 7 embolismic years in each 19-year cycle
    assert_eq!(leap, 7 * cycles);
    // mean year of about 365.2468 days
    let mean = total_days as f64 / (19 * cycles) as f64;
    assert!((mean - 365.2468).abs() < 0.01, "{mean}");
    for a in [3, 6, 8, 11, 14, 17, 19, 5700, 5703] {
        assert!(calendar::is_jewish_leap_year(a), "{a}");
    }
    for a in [1, 2, 4, 5, 7, 18, 5701, 5751] {
        assert!(!calendar::is_jewish_leap_year(a), "{a}");
    }
}

#[test]
fn jewish_year_out_of_range() {
    assert!(calendar::jewish_new_year(i16::MIN).is_none());
    assert!(calendar::jewish_year(i16::MIN).is_none());
    assert!(calendar::jewish_year(i16::MAX).is_none());
    assert!(calendar::jewish_year(i16::MAX - 1).is_some());
    assert!(calendar::jewish_new_year(i16::MAX).is_some());
}

// Meeus, example 9.b
#[test]
fn example_9b_islamic_new_year_1421() {
    let d = calendar::date_from_islamic(&IslamicDate {
        year: 1421,
        month: 1,
        day: 1,
    })
    .unwrap();
    assert_eq!(d, greg(2000, Month::Apr, 6.0));
    assert!(!calendar::is_islamic_leap_year(1421));
}

// Meeus, example 9.c
#[test]
fn example_9c_islamic_from_gregorian() {
    let g = greg(1991, Month::Aug, 13.0);
    // the corresponding Julian date
    let j = Date::from_julian_day_in(time::julian_day(&g), CalType::Julian).unwrap();
    assert_eq!(j, jul(1991, Month::July, 31.0));
    let expected = IslamicDate {
        year: 1412,
        month: 2,
        day: 2,
    };
    assert_eq!(calendar::islamic_from_date(&g), Some(expected));
    assert_eq!(calendar::islamic_from_date(&j), Some(expected));
    // the time of day does not matter
    let g = greg(1991, Month::Aug, 13.99);
    assert_eq!(calendar::islamic_from_date(&g), Some(expected));
}

#[test]
fn islamic_epoch() {
    let epoch = jul(622, Month::July, 16.0);
    assert_eq!(time::julian_day(&epoch), calendar::ISLAMIC_EPOCH_JD);
    let first = IslamicDate {
        year: 1,
        month: 1,
        day: 1,
    };
    assert_eq!(calendar::islamic_from_date(&epoch), Some(first));
    assert_eq!(calendar::date_from_islamic(&first), Some(epoch));
    assert_eq!(
        calendar::islamic_from_date(&jul(622, Month::July, 15.0)),
        None
    );
    assert_eq!(
        calendar::islamic_from_date(&jul(-1000, Month::Jan, 1.0)),
        None
    );
}

#[test]
fn islamic_leap_years() {
    let leap = [2, 5, 7, 10, 13, 16, 18, 21, 24, 26, 29];
    for y in 1..=30 {
        assert_eq!(calendar::is_islamic_leap_year(y), leap.contains(&y), "{y}");
        assert_eq!(
            calendar::is_islamic_leap_year(y + 1410),
            leap.contains(&y),
            "{y}"
        );
    }
}

#[test]
fn islamic_invalid_dates() {
    let bad = [
        (0, 1, 1),
        (-5, 1, 1),
        (1421, 0, 1),
        (1421, 13, 1),
        (1421, 1, 0),
        (1421, 1, 31),
        (1421, 2, 30),
        (1421, 12, 30), // 1421 is a common year
    ];
    for (year, month, day) in bad {
        let d = IslamicDate { year, month, day };
        assert_eq!(calendar::date_from_islamic(&d), None, "{d:?}");
    }
    // the last day of a leap year is valid
    assert!(calendar::is_islamic_leap_year(1420));
    let d = IslamicDate {
        year: 1420,
        month: 12,
        day: 30,
    };
    let civil = calendar::date_from_islamic(&d).unwrap();
    assert_eq!(calendar::islamic_from_date(&civil), Some(d));
    assert_eq!(
        time::days_between(
            &civil,
            &calendar::date_from_islamic(&IslamicDate {
                year: 1421,
                month: 1,
                day: 1
            })
            .unwrap()
        ),
        1.0
    );
}

// Both conversions agree with the arithmetical Islamic calendar of Reingold
// & Dershowitz on every day of 10000 Islamic years. Meeus' own adjustment of
// the day number (JJ) fails on some of these days; see `islamic_from_date`.
#[test]
fn islamic_conversions_match_reference_every_day() {
    let start = calendar::ISLAMIC_EPOCH_JD;
    let days = 3_543_670; // about 10000 Islamic years
    for i in 0..days {
        let jd = start + i as f64;
        let civil = Date::from_julian_day(jd).unwrap();
        let (y, m, d) = islamic_from_rd((jd - 1721424.5) as i64);
        let expected = IslamicDate {
            year: y as i16,
            month: m as u8,
            day: d as u8,
        };
        assert_eq!(
            calendar::islamic_from_date(&civil),
            Some(expected),
            "{civil:?}"
        );
        if d == 1 || d >= 29 {
            assert_eq!(
                calendar::date_from_islamic(&expected),
                Some(civil),
                "{expected:?}"
            );
        }
    }
}

// Days where Meeus' original adjustment goes wrong (found by comparison with
// the reference calendar).
#[test]
fn islamic_from_date_book_adjustment_failures() {
    let cases = [
        (jul(672, Month::Jan, 1.0), (51, 12, 24)),
        (jul(705, Month::Jan, 1.0), (85, 12, 29)),
        (jul(770, Month::Dec, 24.0), (154, 1, 1)),
        (jul(1063, Month::Dec, 25.0), (456, 1, 1)),
    ];
    for (date, (year, month, day)) in cases {
        assert_eq!(
            calendar::islamic_from_date(&date),
            Some(IslamicDate { year, month, day }),
            "{date:?}"
        );
    }
}

#[test]
fn islamic_round_trip_and_range() {
    for year in [1_i16, 2, 29, 30, 31, 1000, 1412, 1446, 5000, 30000] {
        for month in 1..=12 {
            for day in [1, 15, 29] {
                let d = IslamicDate { year, month, day };
                let civil = calendar::date_from_islamic(&d).unwrap();
                assert_eq!(calendar::islamic_from_date(&civil), Some(d));
            }
        }
    }
    // the civil year of a very late Islamic year still fits; the reverse
    // direction reports years beyond i16
    assert!(
        calendar::date_from_islamic(&IslamicDate {
            year: i16::MAX,
            month: 12,
            day: 29
        })
        .is_some()
    );
    assert_eq!(
        calendar::islamic_from_date(&greg(i16::MAX, Month::Dec, 31.0)),
        None
    );
}

// Review R1: the docs used to say that 1 Tishri always falls in September
// or October; that is not so in many Julian years, nor in the far future.
#[test]
fn jewish_new_year_month_range() {
    let mut august = 0;
    for a in 3762_i16..=11800 {
        let d = calendar::jewish_new_year(a).unwrap();
        assert_eq!(d.year as i32, a as i32 - 3761);
        match d.month {
            Month::Aug => {
                assert_eq!(d.cal_type, CalType::Julian, "{d:?}");
                assert!(d.decimal_day >= 25.0, "{d:?}");
                august += 1;
            }
            Month::Sept | Month::Oct => {}
            _ => assert!(d.year > 8100, "{d:?}"),
        }
    }
    assert!(august > 200);
    assert_eq!(
        calendar::jewish_new_year(5261).unwrap(),
        jul(1500, Month::Aug, 25.0)
    );
}
