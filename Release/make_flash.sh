make clean
make all
#arm-none-eabi-objcopy -O binary stm32f401ccAudioNative.elf stm32f401ccAudioNative.bin
#st-flash write stm32f401ccAudioNative.bin 0x8000000
#openocd -f interface/stlink.cfg -f target/stm32f4x.cfg -c "program stm32f401ccAudioNative.bin exit 0x08000000"
openocd -f interface/stlink.cfg -f target/stm32f4x.cfg -c "program stm32f401ccAudioNative.elf verify reset exit"

arm-none-eabi-gcc --version
openocd --version

