use crate::tables::{self, Envelope, DC_FILTER_SIZE, DECIMATE_FACTOR, FIR_SIZE, TONE_CHANNELS};

#[derive(Debug, Clone, Copy, Default)]
pub struct ToneChannel {
    tone_period: u32,
    tone_counter: u32,
    tone: u32,
    /// Whether the tone envelope is off.
    t_off: bool,
    /// Whether the noise envelope is off.
    n_off: bool,
    e_on: bool,
    volume: u32,
    pan_left: f64,
    pan_right: f64,
}

impl ToneChannel {
    pub fn set_tone(&mut self, mut period: u32) {
        period &= 0xfff;
        self.tone_period = (period == 0) as u32 | period;
    }

    pub fn set_mixer(&mut self, t_off: bool, n_off: bool, e_on: bool) {
        self.t_off = t_off;
        self.n_off = n_off;
        self.e_on = e_on;
    }

    pub fn set_pan(&mut self, pan: f64, is_eqp: bool) {
        if is_eqp {
            self.pan_left = (1.0 - pan).sqrt();
            self.pan_right = (pan).sqrt();
        } else {
            self.pan_left = 1.0 - pan;
            self.pan_right = pan;
        }
    }

    pub fn set_volume(&mut self, volume: u32) {
        self.volume = volume & 0xf;
    }

    pub fn update_tone(&mut self) -> u32 {
        self.tone_counter += 1;
        if self.tone_counter >= self.tone_period {
            self.tone_counter = 0;
            self.tone ^= 1;
        }
        self.tone
    }
}

#[derive(Debug, Clone, Copy, Default)]
struct Interpolator {
    c: [f64; 4],
    y: [f64; 4],
}

#[derive(Debug, Clone, Copy)]
struct DcFilter {
    sum: f64,
    delay: [f64; DC_FILTER_SIZE],
}

impl Default for DcFilter {
    fn default() -> Self {
        Self {
            sum: 0.0,
            delay: [0.0; DC_FILTER_SIZE],
        }
    }
}

impl DcFilter {
    fn dc_filter(&mut self, index: usize, x: f64) -> f64 {
        self.sum += -self.delay[index] + x;
        self.delay[index] = x;
        x - self.sum / DC_FILTER_SIZE as f64
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Ayumi {
    channels: [ToneChannel; TONE_CHANNELS],
    noise_period: u32,
    noise_counter: u32,
    noise: u32,
    envelope_counter: u32,
    envelope_period: u32,
    envelope_shape: u32,
    envelope_segment: u32,
    envelope: u32,
    step: f64,
    x: f64,
    is_ym: bool,
    left: f64,
    right: f64,
    fir_left: [f64; FIR_SIZE * 2],
    fir_right: [f64; FIR_SIZE * 2],
    fir_index: usize,
    interpolator_left: Interpolator,
    interpolator_right: Interpolator,
    dc_left: DcFilter,
    dc_right: DcFilter,
    dc_index: usize,
}

impl Default for Ayumi {
    fn default() -> Self {
        Self {
            channels: [ToneChannel::default(); TONE_CHANNELS],
            noise_period: 0,
            noise_counter: 0,
            noise: 0,
            envelope_counter: 0,
            envelope_period: 0,
            envelope_shape: 0,
            envelope_segment: 0,
            envelope: 0,
            step: 0.0,
            x: 0.0,
            is_ym: false,
            left: 0.0,
            right: 0.0,
            fir_left: [0.0; FIR_SIZE * 2],
            fir_right: [0.0; FIR_SIZE * 2],
            fir_index: 0,
            interpolator_left: Interpolator::default(),
            interpolator_right: Interpolator::default(),
            dc_left: DcFilter::default(),
            dc_right: DcFilter::default(),
            dc_index: 0,
        }
    }
}

impl Ayumi {
    pub fn new(is_ym: bool, clock_rate: f64, sr: u32) -> Self {
        let mut ayumi = Self {
            step: clock_rate / (sr * 8 * DECIMATE_FACTOR) as f64,
            noise: 1,
            is_ym,
            ..Default::default()
        };
        ayumi.set_envelope(1);
        for chn in &mut ayumi.channels {
            chn.set_tone(1);
        }
        ayumi
    }

    pub fn set_envelope(&mut self, mut period: u32) {
        period &= 0xffff;
        self.envelope_period = (period == 0) as u32 | period;
    }

    pub fn set_envelope_shape(&mut self, shape: u32) {
        self.envelope_shape = shape & 0xf;
        self.envelope_counter = 0;
        self.envelope_segment = 0;
        self.reset_segment();
    }

    pub fn set_noise(&mut self, mut period: u32) {
        period &= 0x1f;
        self.noise_period = (period == 0) as u32 | period;
    }

    pub fn update_noise(&mut self) -> u32 {
        self.noise_counter += 1;
        if self.noise_counter >= (self.noise_period << 1) {
            self.noise_counter = 0;
            let bit0x3 = (self.noise ^ (self.noise >> 3)) & 1;
            self.noise = (self.noise >> 1) | (bit0x3 << 16);
        }
        self.noise & 1
    }

    pub fn slide_up(&mut self) {
        self.envelope += 1;
        if self.envelope > 31 {
            self.envelope_segment ^= 1;
            self.reset_segment()
        }
    }

    pub fn slide_down(&mut self) {
        if self.envelope == 0 {
            self.envelope_segment ^= 1;
            self.reset_segment();
        } else {
            self.envelope -= 1;
        }
    }

    #[inline]
    const fn get_envelope(&self) -> Envelope {
        tables::ENVELOPES[self.envelope_shape as usize][self.envelope_segment as usize]
    }

    pub fn reset_segment(&mut self) {
        if matches!(self.get_envelope(), Envelope::SlideDown | Envelope::HoldTop) {
            self.envelope = 31;
        } else {
            self.envelope = 0;
        }
    }

    pub fn update_envelope(&mut self) -> u32 {
        self.envelope_counter += 1;
        if self.envelope_counter >= self.envelope_period {
            self.envelope_counter = 0;

            // update the current envelope
            match self.get_envelope() {
                Envelope::SlideUp => self.slide_up(),
                Envelope::SlideDown => self.slide_down(),
                Envelope::HoldTop | Envelope::HoldBottom => (),
            }
        }
        self.envelope
    }

    pub fn update_mixer(&mut self) {
        let (noise, envelope) = (self.update_noise(), self.update_envelope());
        (self.left, self.right) = (0.0, 0.0);
        for chn in &mut self.channels {
            let mut out = (chn.update_tone() | chn.t_off as u32) & (noise | chn.n_off as u32);
            out *= if chn.e_on {
                envelope
            } else {
                chn.volume * 2 + 1
            };
            if self.is_ym {
                self.left += tables::YM_DAC_TABLE[out as usize] * chn.pan_left;
                self.right += tables::YM_DAC_TABLE[out as usize] * chn.pan_right;
            } else {
                self.left += tables::AY_DAC_TABLE[out as usize] * chn.pan_left;
                self.right += tables::AY_DAC_TABLE[out as usize] * chn.pan_right;
            }
        }
    }

    pub fn skip(&mut self) {
        for _ in (DECIMATE_FACTOR - 1)..=0 {
            self.x += self.step;
            if self.x >= 1.0 {
                self.x -= 1.0;
                self.update_mixer();
            }
        }
    }

    pub fn process(&mut self) {
        let c_left = self.interpolator_left.c.clone();
        let y_left = self.interpolator_left.y.clone();
        let c_right = self.interpolator_right.c.clone();
        let y_right = self.interpolator_right.y.clone();
        self.fir_index = (self.fir_index + 1) % (FIR_SIZE / DECIMATE_FACTOR as usize - 1);
        for i in DECIMATE_FACTOR as usize - 1..=0 {
            self.x += self.step;
            if self.x >= 1.0 {
                self.x -= 1.0;
                self.interpolator_left.y[0] = y_left[1];
                self.interpolator_left.y[1] = y_left[2];
                self.interpolator_left.y[2] = y_left[3];
                self.interpolator_right.y[0] = y_right[1];
                self.interpolator_right.y[1] = y_right[2];
                self.interpolator_right.y[2] = y_right[3];
                self.update_mixer();
                self.interpolator_left.y[3] = self.left;
                self.interpolator_left.y[3] = self.right;
                let y1 = y_left[2] - y_left[0];
                self.interpolator_left.c[0] = 0.5 * y_left[1] + 0.25 * (y_left[0] + y_left[2]);
                self.interpolator_left.c[1] = 0.5 * y1;
                self.interpolator_left.c[2] = 0.25 * (y_left[3] - y_left[1] - y1);
                let y1 = y_right[2] - y_right[0];
                self.interpolator_right.c[0] = 0.5 * y_right[1] + 0.25 * (y_right[0] + y_right[2]);
                self.interpolator_right.c[1] = 0.5 * y1;
                self.interpolator_right.c[2] = 0.25 * (y_right[3] - y_right[1] - y1);
            }

            let fir_left =
                &mut self.fir_left[FIR_SIZE - self.fir_index * DECIMATE_FACTOR as usize..];
            let fir_right =
                &mut self.fir_right[FIR_SIZE - self.fir_index * DECIMATE_FACTOR as usize..];

            fir_left[i] = (c_left[2] * self.x + c_left[1]) * self.x + c_left[0];
            fir_right[i] = (c_right[2] * self.x + c_right[1]) * self.x + c_right[0];
        }

        let fir_left = &mut self.fir_left[FIR_SIZE - self.fir_index * DECIMATE_FACTOR as usize..];
        let fir_right = &mut self.fir_right[FIR_SIZE - self.fir_index * DECIMATE_FACTOR as usize..];

        self.left = tables::decimate(fir_left);
        self.right = tables::decimate(fir_right);
    }

    pub fn remove_dc(&mut self) {
        self.left = self.dc_left.dc_filter(self.dc_index, self.left);
        self.right = self.dc_right.dc_filter(self.dc_index, self.right);
        self.dc_index = (self.dc_index + 1) & (DC_FILTER_SIZE - 1);
    }
}
