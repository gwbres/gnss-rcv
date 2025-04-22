use gnss_rtk::prelude::{Epoch, Time, TimeOffset};

// Implement if you ever have more external information.
// Starting from v0.8, gnss_rtk offers possibility to provide measurements
// in any supported timescale, and express temporal solution in any supported timescale,
// provided some external information.
// For example, providing |GST-GPST| and |GST-UTC| allows:
//  - precise mixed GPS+GAL navigation,
//  - providing observations in any of GST, GPST or UTC
//  - solving in either GST, GPST or UTC
pub struct NullTime {}

impl Time for NullTime {
    // Provide |BDT-GPST| update
    fn bdt_gpst_time_offset(&mut self, _now: Epoch) -> Option<TimeOffset> {
        None
    }

    // Provide |BDT-GST| update
    fn bdt_gst_time_offset(&mut self, _now: Epoch) -> Option<TimeOffset> {
        None
    }

    // Provide |BDT-UTC| update
    fn bdt_utc_time_offset(&mut self, _now: Epoch) -> Option<TimeOffset> {
        None
    }

    // Provide |GPST-UTC| update
    fn gpst_utc_time_offset(&mut self, _now: Epoch) -> Option<TimeOffset> {
        None
    }

    // Provide |GST-GPST| update
    fn gst_gpst_time_offset(&mut self, _now: Epoch) -> Option<TimeOffset> {
        None
    }

    // Provide |GST-UTC| update
    fn gst_utc_time_offset(&mut self, _now: Epoch) -> Option<TimeOffset> {
        None
    }
}
