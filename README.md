# stm32f401ccAudioNative
stm32 black pill usb sound card

https://www.youtube.com/watch?v=0MmWp3HdV2A

https://www.youtube.com/watch?v=GbiTxVYopDI


High quality, low noise  DAC based on 2 PWM timer channels with virtual software Sigma Delta ADC between stream from usb and PWM output.
There is implemented "sigma-delta floating point encoder" workaround of native stm32f401 limit 10.5 bits on 44100 Hz (1904 levels=84MHz/44.1KHz )
So, we can have for only  $3 ,very low noise , high sound quality solution, which better then most onboard sound cards !

Same technology can be used for esp32 high quality sound rendering.



