use gnss_rtk::prelude::{
    Bias,
    BiasRuntime,
    //IonosphereModel,
    //KbModel,
    TroposphereModel,
};

pub struct BiasSource {}

impl Bias for BiasSource {
    fn ionosphere_bias_m(&self, _: &BiasRuntime) -> f64 {
        // Ionosphere bias (meters)
        //  - build desired model
        //  - or implement completely
        //  - or provide data directly
        0.0
        // IonosphereModel::KbModel(KbModel {
        //     alpha: (0.0, 0.0, 0.0, 0.0),
        //     beta: (0.0, 0.0, 0.0, 0.0),
        //     h_km: 0.0,
        // })
        // .bias_m(rtm)
    }

    fn troposphere_bias_m(&self, rtm: &BiasRuntime) -> f64 {
        // Troposphere bias (meters)
        //  - select desired model
        //  - or implement one yourself
        //  - or provide data directly
        TroposphereModel::Niel.bias_m(rtm)
    }
}
