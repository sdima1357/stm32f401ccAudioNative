openocd -f interface/stlink.cfg -f target/stm32f4x.cfg -c "program stm32f401ccAudioNative.elf verify reset exit"
