/* USER CODE BEGIN Header */
/**
  ******************************************************************************
  * @file           : main.h
  * @brief          : Header for main.c file.
  *                   This file contains the common defines of the application.
  ******************************************************************************
  * @attention
  *
  * <h2><center>&copy; Copyright (c) 2022 STMicroelectronics.
  * All rights reserved.</center></h2>
  *
  * This software component is licensed by ST under Ultimate Liberty license
  * SLA0044, the "License"; You may not use this file except in compliance with
  * the License. You may obtain a copy of the License at:
  *                             www.st.com/SLA0044
  *
  ******************************************************************************
  */
/* USER CODE END Header */

/* Define to prevent recursive inclusion -------------------------------------*/
#ifndef __MAIN_H
#define __MAIN_H

#ifdef __cplusplus
extern "C" {
#endif

/* Includes ------------------------------------------------------------------*/
#include "stm32f4xx_hal.h"

#include "stm32f4xx_ll_tim.h"
#include "stm32f4xx_ll_bus.h"
#include "stm32f4xx_ll_cortex.h"
#include "stm32f4xx_ll_rcc.h"
#include "stm32f4xx_ll_system.h"
#include "stm32f4xx_ll_utils.h"
#include "stm32f4xx_ll_pwr.h"
#include "stm32f4xx_ll_gpio.h"
#include "stm32f4xx_ll_dma.h"

#include "stm32f4xx_ll_exti.h"

/* Private includes ----------------------------------------------------------*/
/* USER CODE BEGIN Includes */

/* USER CODE END Includes */

/* Exported types ------------------------------------------------------------*/
/* USER CODE BEGIN ET */

/* USER CODE END ET */

/* Exported constants --------------------------------------------------------*/
/* USER CODE BEGIN EC */

/* USER CODE END EC */

/* Exported macro ------------------------------------------------------------*/
/* USER CODE BEGIN EM */

/* USER CODE END EM */

/* Exported functions prototypes ---------------------------------------------*/
void Error_Handler(void);

/* USER CODE BEGIN EFP */

/* USER CODE END EFP */

/* Private defines -----------------------------------------------------------*/
#define LED_Pin GPIO_PIN_13
#define LED_GPIO_Port GPIOC
#define KEY_Pin GPIO_PIN_0
#define KEY_GPIO_Port GPIOA
#define LCD1_CMD_Pin GPIO_PIN_4
#define LCD1_CMD_GPIO_Port GPIOA
#define MUTE_OUT_Pin GPIO_PIN_2
#define MUTE_OUT_GPIO_Port GPIOB
#define LCD_CMD_Pin GPIO_PIN_4
#define LCD_CMD_GPIO_Port GPIOB
#define LCD_RESET_Pin GPIO_PIN_8
#define LCD_RESET_GPIO_Port GPIOB
#define LCD_BL_Pin GPIO_PIN_9
#define LCD_BL_GPIO_Port GPIOB
/* USER CODE BEGIN Private defines */
/*
//!!!!!!!!!!!!!stm32f401ccu6!!!!!!!!!!!!!
V5+			VIN
GND			GND
B12   (i2s2_ws)         LCK
B15   (i2s2_sd)         DIN
B10   (i2s2_ck)         BCK
GND                     SCK
PB2   (mute_out)        XSMT

!!!V3.3+,GND ,LCD_BL  & LCD_RESET are common for both displays!!!

LEFT or siigle display:
GND                     GND
V3.3+                   VCC
PB3 (SPI3_SCK)      	SCK            
PB5 (SPI3_MOSI)    	SDA
PB8 (LCD_RESET)  	RES 
PB4 (LCD_CMD)     	DC     (data/command)
PB9(LCD_BL)          	BLK  (back light)  

RIGHT  display(optional):
GND                     GND
V3.3+                   VCC
PA5 (SPI1_SCK)      	SCK            
PA7 (SPI1_MOSI)    	SDA
PB8 (LCD_RESET)  	RES 
PA4 (LCD1_CMD)     	DC     (data/command)
PB9(LCD_BL)          	BLK  (back light)  

rotate encoder (optional)   
GND                     GND
V3.3+                   +
PA0                     SW
PA15                    DT
PA1                     CLK


*/
/* USER CODE END Private defines */

#ifdef __cplusplus
}
#endif

#endif /* __MAIN_H */

/************************ (C) COPYRIGHT STMicroelectronics *****END OF FILE****/
