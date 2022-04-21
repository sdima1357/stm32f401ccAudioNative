/* USER CODE BEGIN Header */
/**
  ******************************************************************************
  * @file           : main.c
  * @brief          : Main program body
  ******************************************************************************
  * @attention
  *
  * <h2><center>&copy; Copyright (c) 2022 STMicroelectronics.
  * All rights reserved.</center></h2>
  *
  * This software component is licensed by ST under BSD 3-Clause license,
  * the "License"; You may not use this file except in compliance with the
  * License. You may obtain a copy of the License at:
  *                        opensource.org/licenses/BSD-3-Clause
  *
  ******************************************************************************
  */
/* USER CODE END Header */
/* Includes ------------------------------------------------------------------*/
#include "main.h"
#include "usb_device.h"

/* Private includes ----------------------------------------------------------*/
/* USER CODE BEGIN Includes */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "st7789.h"
unsigned char *getImageData();
/* USER CODE END Includes */

/* Private typedef -----------------------------------------------------------*/
/* USER CODE BEGIN PTD */


/* USER CODE END PTD */

/* Private define ------------------------------------------------------------*/
/* USER CODE BEGIN PD */

//**********************start user define************!!!!!
// cpu clock 84 for stm32f401
#define SYS_CLK_MHZ 84
// use internal timers for sound out
#define USE_PWM

// from 0 to 2 SPI LCD
#define NUMBER_OF_LCD 1

//slower for float - faster for integer
//#define USE_FLOAT_SIGMA

// integration order 1 or 2
//#define ORDER 1
#define ORDER 1


// if disabled - direct output (order 0)
#define SIGMA_DELTA

// if enabled - use bilinear interpolation
#define BILINEAR

//**********************end user define area************!!!!!!!!


#ifdef  USE_PWM
#define MAX_LEVELS (SYS_CLK_MHZ*1000000/44100)
#define MAX_VOL (((int)(MAX_LEVELS/SFACT))&~1)
#ifdef  SIGMA_DELTA
#define SFACT 3.1
#else
#define SFACT 1.0
#endif
#else   //No PWM, I2S output to external DAC
#define MAX_VOL (65536)
#define SFACT 1.0
#endif

#define SIGMA_BITS  16
#define SIGMA       (1<<SIGMA_BITS)


#define USB_DATA_BITS   16
#define USB_DATA_BITS_H (USB_DATA_BITS-1)

#define N_SIZE_BITS (8)
#define N_SIZE (1<<N_SIZE_BITS)

#define ASBUF_SIZE     (N_SIZE*2)
int16_t baudio_buffer[ASBUF_SIZE];


#define TIME_BIT_SCALE_FACT 18u
#define TIME_SCALE_FACT     (1u<<TIME_BIT_SCALE_FACT)

// bilinear is more precision
#ifdef BILINEAR
#define BIT_SHIFT_SCALE_FACT 12
#define SHIFT_SCALE_FACT  (1<<BIT_SHIFT_SCALE_FACT)
#endif


//#define MAX_VOL (2046/TIME_SCALE_FACT)


/* USER CODE END PD */

/* Private macro -------------------------------------------------------------*/
/* USER CODE BEGIN PM */

/* USER CODE END PM */

/* Private variables ---------------------------------------------------------*/
I2S_HandleTypeDef hi2s2;
DMA_HandleTypeDef hdma_spi2_tx;

SPI_HandleTypeDef hspi1;
SPI_HandleTypeDef hspi3;

TIM_HandleTypeDef htim3;

UART_HandleTypeDef huart1;

/* USER CODE BEGIN PV */

/* USER CODE END PV */

/* Private function prototypes -----------------------------------------------*/
void SystemClock_Config(void);
static void MX_GPIO_Init(void);
static void MX_DMA_Init(void);
static void MX_I2S2_Init(void);
static void MX_SPI1_Init(void);
static void MX_SPI3_Init(void);
static void MX_TIM1_Init(void);
static void MX_USART1_UART_Init(void);
static void MX_TIM3_Init(void);
/* USER CODE BEGIN PFP */

/* USER CODE END PFP */

/* Private user code ---------------------------------------------------------*/
/* USER CODE BEGIN 0 */

uint16_t* VoiceBuff0 = (uint16_t*)&baudio_buffer[0];
uint16_t* VoiceBuff1 = (uint16_t*)&baudio_buffer[N_SIZE];


void TIM1_TE1()
{
}
void TIM1_TC2()
{
}
void TIM1_HT2()
{
}
void TIM1_TE2()
{
}


void init_timers()
{

	  TIM1->ARR = MAX_VOL+1;
	  LL_AHB1_GRP1_EnableClock(LL_AHB1_GRP1_PERIPH_DMA2);
	  LL_DMA_ConfigAddresses(DMA2, LL_DMA_STREAM_1, (uint32_t)&VoiceBuff0[0], (uint32_t)&TIM1->CCR1, LL_DMA_GetDataTransferDirection(DMA2, LL_DMA_STREAM_1));
	  LL_DMA_SetDataLength(DMA2, LL_DMA_STREAM_1, N_SIZE);
	  LL_DMA_EnableIT_TC(DMA2, LL_DMA_STREAM_1);
	  LL_DMA_EnableIT_TE(DMA2, LL_DMA_STREAM_1);
	  LL_DMA_EnableIT_HT(DMA2, LL_DMA_STREAM_1);

	  LL_DMA_ConfigAddresses(DMA2, LL_DMA_STREAM_2, (uint32_t)&VoiceBuff1[0], (uint32_t)&TIM1->CCR2, LL_DMA_GetDataTransferDirection(DMA2, LL_DMA_STREAM_2));
	  LL_DMA_SetDataLength(DMA2, LL_DMA_STREAM_2, N_SIZE);
	  LL_DMA_EnableIT_TC(DMA2, LL_DMA_STREAM_2);
	  LL_DMA_EnableIT_HT(DMA2, LL_DMA_STREAM_2);
	  LL_DMA_EnableIT_TE(DMA2, LL_DMA_STREAM_2);



	  /***************************/
	  /* Enable the DMA transfer */
	  /***************************/
	  LL_DMA_EnableStream(DMA2, LL_DMA_STREAM_1);
	  LL_DMA_EnableStream(DMA2, LL_DMA_STREAM_2);

	  LL_TIM_EnableDMAReq_CC1(TIM1);
	  LL_TIM_EnableDMAReq_CC2(TIM1);

	  LL_TIM_CC_EnableChannel(TIM1, LL_TIM_CHANNEL_CH1);
	  LL_TIM_CC_EnableChannel(TIM1, LL_TIM_CHANNEL_CH2);



	  LL_TIM_OC_EnablePreload(TIM1, LL_TIM_CHANNEL_CH1);
	  LL_TIM_OC_EnablePreload(TIM1, LL_TIM_CHANNEL_CH2);

	  LL_TIM_CC_EnableChannel(TIM1, LL_TIM_CHANNEL_CH1N);
	  LL_TIM_CC_EnableChannel(TIM1, LL_TIM_CHANNEL_CH2N);
#if 0
	  //LL_TIM_OC_SetPolarity(TIM1, LL_TIM_CHANNEL_CH1,LL_TIM_OCPOLARITY_HIGH);
	  //LL_TIM_OC_SetPolarity(TIM1, LL_TIM_CHANNEL_CH1N,LL_TIM_OCPOLARITY_LOW);
	  LL_TIM_OC_SetPolarity(TIM1, LL_TIM_CHANNEL_CH1,LL_TIM_OCPOLARITY_HIGH);
	  LL_TIM_OC_SetPolarity(TIM1, LL_TIM_CHANNEL_CH1N,LL_TIM_OCPOLARITY_HIGH);
	  LL_TIM_OC_SetPolarity(TIM1, LL_TIM_CHANNEL_CH2,LL_TIM_OCPOLARITY_HIGH);
	  LL_TIM_OC_SetPolarity(TIM1, LL_TIM_CHANNEL_CH2N,LL_TIM_OCPOLARITY_HIGH);
      int res;
	  res = LL_TIM_OC_GetPolarity(TIM1, LL_TIM_CHANNEL_CH1);
	  printf("OC_GetPolarity %d\n",res);
	  res = LL_TIM_OC_GetPolarity(TIM1, LL_TIM_CHANNEL_CH1N);
	  printf("OC_GetPolarity %d\n",res);
#endif


	  LL_TIM_EnableCounter(TIM1);

	  LL_TIM_GenerateEvent_UPDATE(TIM1);
	  LL_TIM_GenerateEvent_CC1(TIM1);
	  

}
//void init_timer3LL()
//{
//	LL_TIM_EnableCounter(TIM3);
//	LL_TIM_GenerateEvent_UPDATE(TIM3);
//	LL_TIM_GenerateEvent_CC1(TIM3);
//}

uint32_t AUDIO_PeriodicTC_FS_Counter;
void retarget_put_char(char p)
{
	HAL_UART_Transmit(&huart1, &p, 1, 0xffffff); // send message via UART
}
int _write(int fd, char* ptr, int len)
{
    (void)fd;
    int i = 0;
    while (ptr[i] && (i < len))
    {
    	if (ptr[i] == '\r')
    	{

    	}
    	else
    	{
			retarget_put_char((int)ptr[i]);
			if (ptr[i] == '\n')
			{
				retarget_put_char((int)'\r');
			}
    	}
        i++;
    }
    return len;
}







volatile int AUDIO_OUT_Play_Counter;
volatile int AUDIO_OUT_ChangeBuffer_Counter;
volatile int TransferComplete_CallBack_FS_Counter;
volatile int HalfTransfer_CallBack_FS_Counter;



//volatile fl = 0;
uint32_t   samplesInBuff        = 0;
int32_t   samplesInBuffH       = 0;
int16_t * usb_SndBuffer        = 0;


uint32_t  samplesInBuffScaled  = 0;
uint32_t  readPositionXScaled  = 0;  //readPosition<<12;
uint32_t  readSpeedXScaled     = 1.01*TIME_SCALE_FACT*(48000.0/(44100*SFACT));


struct LR
{
	int16_t L;
	int16_t R;
};
volatile int32_t ampL = 0;
volatile int32_t ampR = 0;
volatile int32_t LLL = 0;
volatile int32_t RRR = 0;
volatile int UsbSamplesAvail = 0;
struct LR getNextSampleLR()
{
	struct LR res ;
	res.L = 0;
	res.R = 0;
	//if(usb_SndBuffer)
	HAL_GPIO_WritePin(LED_GPIO_Port,LED_Pin,GPIO_PIN_SET);

	if(UsbSamplesAvail)
	{
		HAL_GPIO_WritePin(LED_GPIO_Port,LED_Pin,GPIO_PIN_RESET);
		readPositionXScaled += readSpeedXScaled;
		if(readPositionXScaled >= samplesInBuffScaled)
		{
			readPositionXScaled -= samplesInBuffScaled;
		}

		uint32_t readPositionIntC = readPositionXScaled>>TIME_BIT_SCALE_FACT;
		uint32_t readPositionIntN = readPositionIntC+1;


		int16_t L  = usb_SndBuffer[readPositionIntC*2+0];//+usb_SndBuffer[readPositionIntN*2+0];
		int16_t R  = usb_SndBuffer[readPositionIntC*2+1];//+usb_SndBuffer[readPositionIntN*2+1];

#ifdef BILINEAR
		if(readPositionIntN>=samplesInBuff)
		{
			readPositionIntN -= samplesInBuff;
		}
		int16_t LS  = usb_SndBuffer[readPositionIntN*2+0];//+usb_SndBuffer[readPositionIntN*2+0];
		int16_t RS  = usb_SndBuffer[readPositionIntN*2+1];//+usb_SndBuffer[readPositionIntN*2+1];

		int32_t W1  =  readPositionXScaled>>(TIME_BIT_SCALE_FACT-BIT_SHIFT_SCALE_FACT) & ((1<<BIT_SHIFT_SCALE_FACT)-1);
		int32_t W  =  (1<<BIT_SHIFT_SCALE_FACT)-W1;
		res.L = (L*W+LS*W1)>>BIT_SHIFT_SCALE_FACT;
		res.R = (R*W+RS*W1)>>BIT_SHIFT_SCALE_FACT;
#else  //floor
		res.L = L;
		res.R = R;
#endif
	}
    return res;
}
void  AUDIO_Init(uint32_t AudioFreq, uint32_t Volume, uint32_t options)
{
	printf("AUDIO_Init %d %d %d\n",AudioFreq,Volume,options);
}

void AUDIO_OUT_Start(uint16_t* pBuffer, uint32_t Size)
{
	printf("AUDIO_OUT_Start %p Size = %d\n",pBuffer,Size);
	samplesInBuff = Size/4; //(L & R channels short)
	samplesInBuffH = samplesInBuff/2;

	usb_SndBuffer = pBuffer;
	AUDIO_OUT_Play_Counter++;

	samplesInBuffScaled = samplesInBuff<<TIME_BIT_SCALE_FACT;
	readPositionXScaled = samplesInBuffScaled/2;

}

uint16_t  lastDmaAccessTime;
int       lastDmaPos;
uint32_t   outDmaSpeedScaled;
int       sForcedSpeed;
uint16_t  lastAudioUsbTimeStamp;

int prevPos = -1 ;

void AUDIO_OUT_Periodic(uint16_t* pBuffer, uint32_t Size)
{
	AUDIO_PeriodicTC_FS_Counter++;
	if(!usb_SndBuffer) return ;
	int cPos  = (pBuffer - (uint16_t*)usb_SndBuffer)/2;

#if (NUMBER_OF_LCD>0)
	for(int s=0;s<Size/4;s++)
	{
		int16_t LL = pBuffer[s*2];
		int16_t RR = pBuffer[s*2+1];
		if(LL>ampL)
		{
			ampL = LL;
		}
		if(RR>ampR)
		{
			ampR = RR;
		}
	}
#endif

	if(cPos==0)
	{

		uint16_t lastAudioUsbTimeStampNew = TIM3->CNT;
		//
		uint16_t timeForRecivedSamples = lastAudioUsbTimeStampNew - lastAudioUsbTimeStamp;

        lastAudioUsbTimeStamp = lastAudioUsbTimeStampNew;

		uint16_t timeFromLastDMA = lastAudioUsbTimeStampNew - lastDmaAccessTime;


		//where i am ?

		int approximateSamplesOutedFromLastDMA  = ((float)timeFromLastDMA)*outDmaSpeedScaled/TIME_SCALE_FACT;

		int appDistance  =  (int)(lastDmaPos + approximateSamplesOutedFromLastDMA )-cPos;


        if(UsbSamplesAvail)
		{
			int dC = appDistance - prevPos;
			while(dC>samplesInBuffH)  dC-=samplesInBuff;
			while(dC<-samplesInBuffH) dC+=samplesInBuff;
			if(dC > samplesInBuffH/8 || dC<-samplesInBuffH/8 ) //seems completely lost sync , force set frequency
			{
		        float  speedFact =  (((float)timeForRecivedSamples) * outDmaSpeedScaled)/TIME_SCALE_FACT;
		        sForcedSpeed  = (int)(((float)samplesInBuff)*TIME_SCALE_FACT / speedFact);
				readSpeedXScaled = sForcedSpeed;
			}
			else
			{
				//ok - only phase tune
				readSpeedXScaled-=dC;
			}
		}
        else
        {
        	readPositionXScaled = samplesInBuffScaled/2;
        }
		prevPos =  appDistance;
		UsbSamplesAvail = samplesInBuff*((int)SFACT)*8;
	}
	return 0;
}

typedef struct
{
	float x;
	float y;
} vec2f;

typedef float mat4f[4] ;
void printMat(mat4f m2)
{
	printf("%1.3f %1.3f\n",m2[0],m2[1]);
	printf("%1.3f %1.3f\n",m2[2],m2[3]);
}

void printVec(vec2f v)
{
	printf("x=%1.3f y=%1.3f\n",v.x,v.y);
}



vec2f applyMat(mat4f m2,vec2f v)
{
	vec2f res;
	res.x = v.x*m2[0]+v.y*m2[1];
	res.y = v.x*m2[2]+v.y*m2[3];
	return res;
}



void makeRotMat(mat4f xx,float alpha)
{
	float cosA = cosf(alpha);
	float sinA = sinf(alpha);
	xx[0]  = cosA; xx[1]  = -sinA;
	xx[2]  = sinA; xx[3]  = cosA;
}
void AUDIO_OUT_ChangeBuffer(uint16_t *pBuffer, uint16_t Size)
{
	printf("%p SizeC = %d\n",pBuffer,Size);
		AUDIO_OUT_ChangeBuffer_Counter++;
		//fl = 1;
		//HAL_I2S_Transmit_DMA(&hi2s2, pBuffer,Size);
}
void cleanAudioBuffer()
{
#ifdef PWM
	for(int k=0;k<N_SIZE;k++)
	{
		VoiceBuff0[k] = MAX_VOL/2;
		VoiceBuff1[k] = MAX_VOL/2;
	}
#else
	for(int k=0;k<N_SIZE;k++)
	{
		VoiceBuff0[k] = 0;
		VoiceBuff1[k] = 0;
	}
#endif
}

struct sigmaDeltaStorage
{
	float integral;
	float y;
};
struct sigmaDeltaStorage2
{
	float integral0;
	float integral1;
	float y;
};


float funcN1(float x,int N)
{
	//float y = 2.0f*floorf(x*0.5f)+1.0f;
	float y = floorf(x+0.5f);
	//saturate
	if(y>N) y = N;
	if(y<-N) y = -N;
	return y;
}


float sigma_delta_slow(struct sigmaDeltaStorage* st,float x)
{
	st->integral+= x - st->y;
	st->y = funcN1(st->integral,(int)(MAX_VOL/2));
	return st->y;
}



struct  sigmaDeltaStorage_SCALED
{
	int integral;
	int y;
};
struct  sigmaDeltaStorage2_SCALED
{
	int integral0;
	int integral1;
	int y;
};

float sigma_delta(struct sigmaDeltaStorage* st,float x)
{
	st->integral+= x - st->y;
	st->y = floorf(st->integral+0.5f);
	if(st->y<0)st->y = 0;
	if(st->y>MAX_VOL)st->y = MAX_VOL;

	return st->y;
}

float sigma_delta2(struct sigmaDeltaStorage2* st,float x)
{
	st->integral0+= x             - st->y;
	st->integral1+= st->integral0 - st->y;
	st->y = floorf(st->integral1+0.5f);
	if(st->y<0)st->y = 0;
	if(st->y>MAX_VOL)st->y = MAX_VOL;

	return st->y;
}


int sigma_delta_SCALED(struct sigmaDeltaStorage_SCALED* st,int x_SCALED)
{
	st->integral+= x_SCALED - st->y*SIGMA; //post diff
	st->y		 =(st->integral+SIGMA/2)/SIGMA;
	if(st->y < 0)
		st->y = 0;
	if(st->y > MAX_VOL)
		st->y = MAX_VOL;
	return st->y;
}

int sigma_delta2_SCALED(struct sigmaDeltaStorage2_SCALED* st,int x_SCALED)
{
	st->integral0+= x_SCALED      - st->y*SIGMA;
	st->integral1+= st->integral0 - st->y*SIGMA;
	st->y		 =(st->integral1+SIGMA/2)/SIGMA;
	if(st->y < 0)
		st->y = 0;
	if(st->y > MAX_VOL)
		st->y = MAX_VOL;
	return st->y;
}

struct sigmaDeltaStorage static_L_channel;
struct sigmaDeltaStorage static_R_channel;

struct sigmaDeltaStorage_SCALED static_L_channel_SCALED;
struct sigmaDeltaStorage_SCALED static_R_channel_SCALED;

struct sigmaDeltaStorage2 static_L_channel2;
struct sigmaDeltaStorage2 static_R_channel2;

struct sigmaDeltaStorage2_SCALED static_L_channel2_SCALED;
struct sigmaDeltaStorage2_SCALED static_R_channel2_SCALED;


void TIM1_TC1()
{
	uint16_t prevTime;
	uint16_t  tfl;
	prevTime = lastDmaAccessTime;
	lastDmaAccessTime = TIM3->CNT;
	lastDmaPos  = readPositionXScaled>>TIME_BIT_SCALE_FACT;
	tfl      = lastDmaAccessTime -prevTime;
	outDmaSpeedScaled =  TIME_SCALE_FACT*N_SIZE/(tfl);
    if(UsbSamplesAvail > N_SIZE/2)
    {
    	UsbSamplesAvail -= N_SIZE/2;
    }
    else
    {
    	UsbSamplesAvail = 0;
    }
	TransferComplete_CallBack_FS_Counter++;
#ifndef SIGMA_DELTA
	for(int k=0;k<N_SIZE/2;k++)
	{
		struct LR tt =getNextSampleLR(/*k+ASBUF_SIZE/4*/);
		VoiceBuff0[k+N_SIZE/2] = MAX_VOL*(tt.L+(1<<USB_DATA_BITS_H))>>USB_DATA_BITS;
		VoiceBuff1[k+N_SIZE/2] = MAX_VOL*(tt.R+(1<<USB_DATA_BITS_H))>>USB_DATA_BITS;
	}
#else
#ifdef USE_FLOAT_SIGMA
#if    (ORDER==1)
	float a_scale = MAX_VOL/65536.0f;
	for(int k=0;k<N_SIZE/2;k++)
	{
		struct LR tt =getNextSampleLR(/*k+ASBUF_SIZE/4*/);
		VoiceBuff0[k+N_SIZE/2] =( int)(sigma_delta2(&static_L_channel2,a_scale*tt.L+ MAX_VOL/2));
		VoiceBuff1[k+N_SIZE/2] =( int)(sigma_delta2(&static_R_channel2,a_scale*tt.R+ MAX_VOL/2));
	}
#else //order == 2
	float a_scale = MAX_VOL/65536.0f;
	for(int k=0;k<N_SIZE/2;k++)
	{
		struct LR tt =getNextSampleLR(/*k+ASBUF_SIZE/4*/);
		VoiceBuff0[k+N_SIZE/2] =( int)(sigma_delta(&static_L_channel,a_scale*tt.L+ MAX_VOL/2));
		VoiceBuff1[k+N_SIZE/2] =( int)(sigma_delta(&static_R_channel,a_scale*tt.R+ MAX_VOL/2));
	}
#endif
#else
#if (ORDER==1)
	for(int k=0;k<N_SIZE/2;k++)
	{
		struct LR tt =getNextSampleLR(/*k+ASBUF_SIZE/4*/);
		VoiceBuff0[k+N_SIZE/2] =sigma_delta_SCALED(&static_L_channel_SCALED,(MAX_VOL*(tt.L+(1<<USB_DATA_BITS_H)))>>(USB_DATA_BITS-SIGMA_BITS));
		VoiceBuff1[k+N_SIZE/2] =sigma_delta_SCALED(&static_R_channel_SCALED,(MAX_VOL*(tt.R+(1<<USB_DATA_BITS_H)))>>(USB_DATA_BITS-SIGMA_BITS));
	}
#else //order == 2
	for(int k=0;k<N_SIZE/2;k++)
	{
		struct LR tt =getNextSampleLR(/*k+ASBUF_SIZE/4*/);
		VoiceBuff0[k+N_SIZE/2] =sigma_delta2_SCALED(&static_L_channel2_SCALED,(MAX_VOL*(tt.L+(1<<USB_DATA_BITS_H)))>>(USB_DATA_BITS-SIGMA_BITS));
		VoiceBuff1[k+N_SIZE/2] =sigma_delta2_SCALED(&static_R_channel2_SCALED,(MAX_VOL*(tt.R+(1<<USB_DATA_BITS_H)))>>(USB_DATA_BITS-SIGMA_BITS));
	}
#endif

#endif
#endif
}
void TIM1_HT1()
{
	HalfTransfer_CallBack_FS_Counter++;
    if(UsbSamplesAvail > N_SIZE/2)
    {
    	UsbSamplesAvail -= N_SIZE/2;
    }
    else
    {
    	UsbSamplesAvail = 0;
    }
	struct LR * bfr = ((struct LR *) baudio_buffer);
#ifndef  SIGMA_DELTA
	for(int k=0;k<N_SIZE/2;k++)
	{
		struct LR tt =getNextSampleLR(/*k+ASBUF_SIZE/4*/);
		VoiceBuff0[k] = MAX_VOL*(tt.L+(1<<USB_DATA_BITS_H))>>USB_DATA_BITS;;
		VoiceBuff1[k] = MAX_VOL*(tt.R+(1<<USB_DATA_BITS_H))>>USB_DATA_BITS;
	}
#else
#ifdef USE_FLOAT_SIGMA
#if (ORDER==1)
	float a_scale = MAX_VOL/65536.0f;
	for(int k=0;k<N_SIZE/2;k++)
	{
		struct LR tt =getNextSampleLR(/*k+ASBUF_SIZE/4*/);
		VoiceBuff0[k] =( int)(sigma_delta2(&static_L_channel2,a_scale*tt.L+ MAX_VOL/2));
		VoiceBuff1[k] =( int)(sigma_delta2(&static_R_channel2,a_scale*tt.R+ MAX_VOL/2));
	}
#else //order == 2
	float a_scale = ((float)MAX_VOL)/(1<<USB_DATA_BITS);
	for(int k=0;k<N_SIZE/2;k++)
	{
		struct LR tt =getNextSampleLR(/*k+ASBUF_SIZE/4*/);
		VoiceBuff0[k] =( int)(sigma_delta(&static_L_channel,a_scale*tt.L+ MAX_VOL/2));
		VoiceBuff1[k] =( int)(sigma_delta(&static_R_channel,a_scale*tt.R+ MAX_VOL/2));
	}
#endif
#else
#if ORDER==1
	for(int k=0;k<N_SIZE/2;k++)
	{
		struct LR tt =getNextSampleLR(/*k+ASBUF_SIZE/4*/);
		VoiceBuff0[k] =sigma_delta_SCALED(&static_L_channel_SCALED,(MAX_VOL*(tt.L+(1<<USB_DATA_BITS_H)))>>(USB_DATA_BITS-SIGMA_BITS));
		VoiceBuff1[k] =sigma_delta_SCALED(&static_R_channel_SCALED,(MAX_VOL*(tt.R+(1<<USB_DATA_BITS_H)))>>(USB_DATA_BITS-SIGMA_BITS));
	}
#else //order == 2
	for(int k=0;k<N_SIZE/2;k++)
	{
		struct LR tt =getNextSampleLR(/*k+ASBUF_SIZE/4*/);
		VoiceBuff0[k] =sigma_delta2_SCALED(&static_L_channel2_SCALED,(MAX_VOL*(tt.L+(1<<USB_DATA_BITS_H)))>>(USB_DATA_BITS-SIGMA_BITS));
		VoiceBuff1[k] =sigma_delta2_SCALED(&static_R_channel2_SCALED,(MAX_VOL*(tt.R+(1<<USB_DATA_BITS_H)))>>(USB_DATA_BITS-SIGMA_BITS));
	}
#endif
#endif
#endif
}
void HAL_I2S_TxCpltCallback(I2S_HandleTypeDef *hi2s)
{
	//fl = 0;
	uint16_t prevTime;
	uint16_t  tfl;
    if(UsbSamplesAvail > N_SIZE/2)
    {
    	UsbSamplesAvail -= N_SIZE/2;
    }
    else
    {
    	UsbSamplesAvail = 0;
    }

	prevTime = lastDmaAccessTime;
	lastDmaAccessTime = TIM3->CNT;
	lastDmaPos  = readPositionXScaled>>TIME_BIT_SCALE_FACT;
	tfl = lastDmaAccessTime -prevTime;
	outDmaSpeedScaled =  TIME_SCALE_FACT*N_SIZE/(tfl);

	TransferComplete_CallBack_FS_Counter++;
    struct LR * bfr = ((struct LR *) baudio_buffer)+ASBUF_SIZE/4;
	for(int k=0;k<ASBUF_SIZE/4;k++)
	{
		bfr[k] = getNextSampleLR(/*k+ASBUF_SIZE/4*/);

	}
}

//void HAL_I2S_TxHalfCpltCallback(I2S_HandleTypeDef *hi2s)
void HAL_I2S_TxHalfCpltCallback(I2S_HandleTypeDef *hi2s)
{
    if(UsbSamplesAvail > N_SIZE/2)
    {
    	UsbSamplesAvail -= N_SIZE/2;
    }
    else
    {
    	UsbSamplesAvail = 0;
    }
	HalfTransfer_CallBack_FS_Counter++;
	struct LR * bfr = ((struct LR *) baudio_buffer);
	for(int k=0;k<ASBUF_SIZE/4;k++)
		bfr[k] = getNextSampleLR(/*k*/);
	//HalfTransfer_CallBack_FS();
}
void HAL_I2S_ErrorCallback(I2S_HandleTypeDef *hi2s)
{
	printf("HAL_I2S_ErrorCallback !!!!!!!!!!!!!!!\n");
  //BSP_AUDIO_OUT_Error_CallBack();
}


#define PI 3.14159265358979323846
#define TAU (2.0 * PI)


/* USER CODE END 0 */

/**
  * @brief  The application entry point.
  * @retval int
  */
int main(void)
{
  /* USER CODE BEGIN 1 */

  /* USER CODE END 1 */

  /* MCU Configuration--------------------------------------------------------*/

  /* Reset of all peripherals, Initializes the Flash interface and the Systick. */
  HAL_Init();

  /* USER CODE BEGIN Init */

  /* USER CODE END Init */

  /* Configure the system clock */
  SystemClock_Config();

  /* USER CODE BEGIN SysInit */

  /* USER CODE END SysInit */

  /* Initialize all configured peripherals */
  MX_GPIO_Init();
  MX_DMA_Init();
  MX_I2S2_Init();
  MX_SPI1_Init();
  MX_SPI3_Init();
  MX_USB_DEVICE_Init();
  MX_TIM1_Init();
  MX_USART1_UART_Init();
  MX_TIM3_Init();
  /* USER CODE BEGIN 2 */

  printf("\nStart program\n");
  cleanAudioBuffer();
  TIM3->PSC  = SYS_CLK_MHZ*10/4-1;
  HAL_TIM_Base_Start(&htim3);
  HAL_StatusTypeDef stat;

#ifdef  USE_PWM
  LL_TIM_DisableAllOutputs(TIM1);
  init_timers();
  HAL_Delay(100);
  LL_TIM_EnableAllOutputs(TIM1);
//  LL_GPIO_InitTypeDef GPIO_InitStruct = {0};
//  GPIO_InitStruct.Pin = LL_GPIO_PIN_8|LL_GPIO_PIN_9;
//  GPIO_InitStruct.Mode = LL_GPIO_MODE_ALTERNATE;
//  GPIO_InitStruct.Speed = LL_GPIO_SPEED_FREQ_LOW;
//  GPIO_InitStruct.OutputType = LL_GPIO_OUTPUT_PUSHPULL;
//  GPIO_InitStruct.Pull = LL_GPIO_PULL_NO;
//  GPIO_InitStruct.Alternate = LL_GPIO_AF_1;
//  LL_GPIO_Init(GPIOA, &GPIO_InitStruct);
#else
  stat = HAL_I2S_Transmit_DMA(&hi2s2, baudio_buffer,ASBUF_SIZE);
  HAL_Delay(100);
#endif
  printf("Start Speed %x\n",readSpeedXScaled);
  printf("sSpeed= %x %04d %04d %04d speed=%x\n",sForcedSpeed,HalfTransfer_CallBack_FS_Counter,TransferComplete_CallBack_FS_Counter,AUDIO_PeriodicTC_FS_Counter,readSpeedXScaled);
#if (NUMBER_OF_LCD)
	  LCD_reset();
	  for(int l=0;l<NUMBER_OF_LCD;l++)
	  {
		  setLCD(l);
		  LCD_init();
		  LCD_setRotation(PORTRAIT_FLIP);
		  {
				uint32_t trt = HAL_GetTick();
				//LCD_fillRect(0, 0, LCD_getWidth(), LCD_getHeight(), rand());
				LCD_fillRectData(0, 0, LCD_getWidth(), LCD_getHeight(),getImageData());
				trt = (HAL_GetTick()-trt);
				printf("SCR_FULL %d uS\n",trt*1000);

		  }
	  }
#endif
  int pcXe_OLD[2]={-1,-1};
  int pcYe_OLD[2]={-1,-1};
  int pcX_OLD[2]={-1,-1};
  int pcY_OLD[2]={-1,-1};
  float elapsed_time = 0;
  int  elapsed_time_ticks = 0;
  int32_t prev_tick = 0;
  int  mode = 0;
  /* USER CODE END 2 */

  /* Infinite loop */
  /* USER CODE BEGIN WHILE */
  int timcnt =0;
  while (1)
  {
    /* USER CODE END WHILE */

    /* USER CODE BEGIN 3 */
	  timcnt++;
	  if(timcnt%100==0)
	  {
	  // printf("%04d %04d %04d %04d %04d\n",HalfTransfer_CallBack_FS_Counter,TransferComplete_CallBack_FS_Counter,AUDIO_PeriodicTC_FS_Counter,AUDIO_OUT_Play_Counter,AUDIO_OUT_ChangeBuffer_Counter);
		  printf("MAX_VOL = %d outSpeed=%0.01f Hz %x elaps = %d mS,ForcedSpeed = %x %04d PLLspeed = %x\n",MAX_VOL,(400000.0f*(float)outDmaSpeedScaled)/TIME_SCALE_FACT,outDmaSpeedScaled,elapsed_time_ticks,sForcedSpeed,AUDIO_PeriodicTC_FS_Counter,readSpeedXScaled);

	  }
	  if(HAL_GPIO_ReadPin(KEY_GPIO_Port,KEY_Pin)==0 && NUMBER_OF_LCD > 0)
	  {
		 if(mode==0)
		 {
			 HAL_GPIO_WritePin(LCD_BL_GPIO_Port, LCD_BL_Pin,GPIO_PIN_RESET);
   		     LCD_reset();
			 HAL_GPIO_WritePin(LCD_BL_GPIO_Port, LCD_BL_Pin,GPIO_PIN_RESET);
   		     HAL_Delay(1000);
			 mode = 1;
		 }
		 else
		 {
			 mode = 0;

			  for(int l=0;l<NUMBER_OF_LCD;l++)
			  {
				  pcXe_OLD[l] = -1;
				  setLCD(l);
				  LCD_init();
				  LCD_setRotation(PORTRAIT_FLIP);
				  LCD_fillRectData(0, 0, LCD_getWidth(), LCD_getHeight(),getImageData());
			  }
			 HAL_GPIO_WritePin(LCD_BL_GPIO_Port, LCD_BL_Pin,GPIO_PIN_SET);
			 HAL_Delay(1000);
		 }
	  }
		//ampL = (ampL*15)/(16);
		//ampR = (ampR*15)/(16);
	  if(mode==0)
	  {
			//LCD_fillRectData(0, 0, LCD_getWidth(), LCD_getHeight(),getImageData());
		    uint32_t trt = HAL_GetTick();
		    int deltaTicks = trt - prev_tick;
		    //300 millisecond
#define decayFact (-1.0f/300)
		    float decay = expf(deltaTicks*decayFact);
		    ampL =ampL*decay;
		    ampR =ampR*decay;
		    prev_tick = trt;

		    for(int l=0;l<NUMBER_OF_LCD;l++)
		    {
		    	setLCD(l);
				vec2f vc = {sqrtf((159-79)*(159-79)+(119-55)*(119-55)),0.0f};
				mat4f mt;
				float ampl = ((l==0)?ampR:ampL)/600.0f/16;
#if (NUMBER_OF_LCD==1)
					ampl = (ampR+ampL)*0.5f/600.0f/16;
#endif
				if(ampl>1.0f) ampl = 1.0f;
				if(ampl<0.0f) ampl = 0.0f;
				makeRotMat(mt,(90+42-84*ampl)/180*((float)M_PI));
				vec2f vcA = applyMat(mt,vc);
				//printMat(mt);
				//printVec(vc);
				///printVec(vcA);
				int pcX = 121;//55
				int pcY = 167;//79
				int pcXe = (int)(pcX+vcA.x);//55
				int pcYe = (int)(pcY-vcA.y);//79
				int pcXs = (int)(pcX+vcA.x*0.22f);
				int pcYs = (int)(pcY-vcA.y*0.22f);
				//if(!(pcXe_OLD[l] == pcXe && pcYe_OLD[l] == pcYe))
				{
					if(pcXe_OLD[l]>0)
					{
						LCD_Draw_LinePNX(pcX_OLD[l]-2-5,pcY_OLD[l],pcXe_OLD[l]-2-5,pcYe_OLD[l],240,4,getImageData());
						LCD_Draw_LinePNX(pcX_OLD[l]-2,pcY_OLD[l],pcXe_OLD[l]-2,pcYe_OLD[l],240,4,getImageData());
					}
					if(pcXe_OLD[l]>0)
					{
					}
					LCD_Draw_LinePNXShadow(pcXs-2-5,pcYs,pcXe-2-5,pcYe,240,4,getImageData());
					//LCD_Draw_LinePNX(pcXs-2-5,pcYs,pcXe-2-5,pcYe,240,4,getImageDarkData());

					uint16_t colors [4];
					colors[0] = color_convertRGB_to16d(0x3e + 20,0x27 + 20,0x2f + 20);
					colors[1] = color_convertRGB_to16d(0xd1/2,0x26/2,0x64/2);
					colors[2] = color_convertRGB_to16d(0xd1,0x26,0x64);
					colors[3] = color_convertRGB_to16d(0x3e,0x27,0x2f);
					LCD_Draw_LineNX(pcXs-2,pcYs,pcXe-2,pcYe,4,(uint8_t*) colors);
				}

//				LCD_Draw_Line(pcX-2,pcY,pcXe-2,pcYe,RED);
//				LCD_Draw_Line(pcX-1,pcY,pcXe-1,pcYe,BLACK);
//				LCD_Draw_Line(pcX,pcY,pcXe,pcYe,BLACK);
//				LCD_Draw_Line(pcX+1,pcY,pcXe+1,pcYe,RED);

				 pcX_OLD[l] = pcXs;
				 pcY_OLD[l] = pcYs;
				 pcXe_OLD[l] = pcXe;
				 pcYe_OLD[l] = pcYe;
		    }
		    elapsed_time_ticks = HAL_GetTick() - trt;
			//printf(txt,"%03d %03d\n",ampL,ampR);
			//LCD_Draw_Text(txt,0,0,GREEN, 2,BLACK);
	  }
	  HAL_Delay(20);

  }
  /* USER CODE END 3 */
}

/**
  * @brief System Clock Configuration
  * @retval None
  */
void SystemClock_Config(void)
{
  RCC_OscInitTypeDef RCC_OscInitStruct = {0};
  RCC_ClkInitTypeDef RCC_ClkInitStruct = {0};

  /** Configure the main internal regulator output voltage
  */
  __HAL_RCC_PWR_CLK_ENABLE();
  __HAL_PWR_VOLTAGESCALING_CONFIG(PWR_REGULATOR_VOLTAGE_SCALE2);
  /** Initializes the RCC Oscillators according to the specified parameters
  * in the RCC_OscInitTypeDef structure.
  */
  RCC_OscInitStruct.OscillatorType = RCC_OSCILLATORTYPE_HSE;
  RCC_OscInitStruct.HSEState = RCC_HSE_ON;
  RCC_OscInitStruct.PLL.PLLState = RCC_PLL_ON;
  RCC_OscInitStruct.PLL.PLLSource = RCC_PLLSOURCE_HSE;
  RCC_OscInitStruct.PLL.PLLM = 25;
  RCC_OscInitStruct.PLL.PLLN = 336;
  RCC_OscInitStruct.PLL.PLLP = RCC_PLLP_DIV4;
  RCC_OscInitStruct.PLL.PLLQ = 7;
  if (HAL_RCC_OscConfig(&RCC_OscInitStruct) != HAL_OK)
  {
    Error_Handler();
  }
  /** Initializes the CPU, AHB and APB buses clocks
  */
  RCC_ClkInitStruct.ClockType = RCC_CLOCKTYPE_HCLK|RCC_CLOCKTYPE_SYSCLK
                              |RCC_CLOCKTYPE_PCLK1|RCC_CLOCKTYPE_PCLK2;
  RCC_ClkInitStruct.SYSCLKSource = RCC_SYSCLKSOURCE_PLLCLK;
  RCC_ClkInitStruct.AHBCLKDivider = RCC_SYSCLK_DIV1;
  RCC_ClkInitStruct.APB1CLKDivider = RCC_HCLK_DIV2;
  RCC_ClkInitStruct.APB2CLKDivider = RCC_HCLK_DIV1;

  if (HAL_RCC_ClockConfig(&RCC_ClkInitStruct, FLASH_LATENCY_2) != HAL_OK)
  {
    Error_Handler();
  }
}

/**
  * @brief I2S2 Initialization Function
  * @param None
  * @retval None
  */
static void MX_I2S2_Init(void)
{

  /* USER CODE BEGIN I2S2_Init 0 */

  /* USER CODE END I2S2_Init 0 */

  /* USER CODE BEGIN I2S2_Init 1 */

  /* USER CODE END I2S2_Init 1 */
  hi2s2.Instance = SPI2;
  hi2s2.Init.Mode = I2S_MODE_MASTER_TX;
  hi2s2.Init.Standard = I2S_STANDARD_PHILIPS;
  hi2s2.Init.DataFormat = I2S_DATAFORMAT_16B;
  hi2s2.Init.MCLKOutput = I2S_MCLKOUTPUT_DISABLE;
  hi2s2.Init.AudioFreq = I2S_AUDIOFREQ_48K;
  hi2s2.Init.CPOL = I2S_CPOL_LOW;
  hi2s2.Init.ClockSource = I2S_CLOCK_PLL;
  hi2s2.Init.FullDuplexMode = I2S_FULLDUPLEXMODE_DISABLE;
  if (HAL_I2S_Init(&hi2s2) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN I2S2_Init 2 */

  /* USER CODE END I2S2_Init 2 */

}

/**
  * @brief SPI1 Initialization Function
  * @param None
  * @retval None
  */
static void MX_SPI1_Init(void)
{

  /* USER CODE BEGIN SPI1_Init 0 */

  /* USER CODE END SPI1_Init 0 */

  /* USER CODE BEGIN SPI1_Init 1 */

  /* USER CODE END SPI1_Init 1 */
  /* SPI1 parameter configuration*/
  hspi1.Instance = SPI1;
  hspi1.Init.Mode = SPI_MODE_MASTER;
  hspi1.Init.Direction = SPI_DIRECTION_2LINES;
  hspi1.Init.DataSize = SPI_DATASIZE_8BIT;
  hspi1.Init.CLKPolarity = SPI_POLARITY_HIGH;
  hspi1.Init.CLKPhase = SPI_PHASE_2EDGE;
  hspi1.Init.NSS = SPI_NSS_SOFT;
  hspi1.Init.BaudRatePrescaler = SPI_BAUDRATEPRESCALER_2;
  hspi1.Init.FirstBit = SPI_FIRSTBIT_MSB;
  hspi1.Init.TIMode = SPI_TIMODE_DISABLE;
  hspi1.Init.CRCCalculation = SPI_CRCCALCULATION_DISABLE;
  hspi1.Init.CRCPolynomial = 10;
  if (HAL_SPI_Init(&hspi1) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN SPI1_Init 2 */

  /* USER CODE END SPI1_Init 2 */

}

/**
  * @brief SPI3 Initialization Function
  * @param None
  * @retval None
  */
static void MX_SPI3_Init(void)
{

  /* USER CODE BEGIN SPI3_Init 0 */

  /* USER CODE END SPI3_Init 0 */

  /* USER CODE BEGIN SPI3_Init 1 */

  /* USER CODE END SPI3_Init 1 */
  /* SPI3 parameter configuration*/
  hspi3.Instance = SPI3;
  hspi3.Init.Mode = SPI_MODE_MASTER;
  hspi3.Init.Direction = SPI_DIRECTION_2LINES;
  hspi3.Init.DataSize = SPI_DATASIZE_8BIT;
  hspi3.Init.CLKPolarity = SPI_POLARITY_HIGH;
  hspi3.Init.CLKPhase = SPI_PHASE_2EDGE;
  hspi3.Init.NSS = SPI_NSS_SOFT;
  hspi3.Init.BaudRatePrescaler = SPI_BAUDRATEPRESCALER_2;
  hspi3.Init.FirstBit = SPI_FIRSTBIT_MSB;
  hspi3.Init.TIMode = SPI_TIMODE_DISABLE;
  hspi3.Init.CRCCalculation = SPI_CRCCALCULATION_DISABLE;
  hspi3.Init.CRCPolynomial = 10;
  if (HAL_SPI_Init(&hspi3) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN SPI3_Init 2 */

  /* USER CODE END SPI3_Init 2 */

}

/**
  * @brief TIM1 Initialization Function
  * @param None
  * @retval None
  */
static void MX_TIM1_Init(void)
{

  /* USER CODE BEGIN TIM1_Init 0 */

  /* USER CODE END TIM1_Init 0 */

  LL_TIM_InitTypeDef TIM_InitStruct = {0};
  LL_TIM_OC_InitTypeDef TIM_OC_InitStruct = {0};
  LL_TIM_BDTR_InitTypeDef TIM_BDTRInitStruct = {0};

  LL_GPIO_InitTypeDef GPIO_InitStruct = {0};

  /* Peripheral clock enable */
  LL_APB2_GRP1_EnableClock(LL_APB2_GRP1_PERIPH_TIM1);

  /* TIM1 DMA Init */

  /* TIM1_CH1 Init */
  LL_DMA_SetChannelSelection(DMA2, LL_DMA_STREAM_1, LL_DMA_CHANNEL_6);

  LL_DMA_SetDataTransferDirection(DMA2, LL_DMA_STREAM_1, LL_DMA_DIRECTION_MEMORY_TO_PERIPH);

  LL_DMA_SetStreamPriorityLevel(DMA2, LL_DMA_STREAM_1, LL_DMA_PRIORITY_HIGH);

  LL_DMA_SetMode(DMA2, LL_DMA_STREAM_1, LL_DMA_MODE_CIRCULAR);

  LL_DMA_SetPeriphIncMode(DMA2, LL_DMA_STREAM_1, LL_DMA_PERIPH_NOINCREMENT);

  LL_DMA_SetMemoryIncMode(DMA2, LL_DMA_STREAM_1, LL_DMA_MEMORY_INCREMENT);

  LL_DMA_SetPeriphSize(DMA2, LL_DMA_STREAM_1, LL_DMA_PDATAALIGN_HALFWORD);

  LL_DMA_SetMemorySize(DMA2, LL_DMA_STREAM_1, LL_DMA_MDATAALIGN_HALFWORD);

  LL_DMA_DisableFifoMode(DMA2, LL_DMA_STREAM_1);

  /* TIM1_CH2 Init */
  LL_DMA_SetChannelSelection(DMA2, LL_DMA_STREAM_2, LL_DMA_CHANNEL_6);

  LL_DMA_SetDataTransferDirection(DMA2, LL_DMA_STREAM_2, LL_DMA_DIRECTION_MEMORY_TO_PERIPH);

  LL_DMA_SetStreamPriorityLevel(DMA2, LL_DMA_STREAM_2, LL_DMA_PRIORITY_HIGH);

  LL_DMA_SetMode(DMA2, LL_DMA_STREAM_2, LL_DMA_MODE_CIRCULAR);

  LL_DMA_SetPeriphIncMode(DMA2, LL_DMA_STREAM_2, LL_DMA_PERIPH_NOINCREMENT);

  LL_DMA_SetMemoryIncMode(DMA2, LL_DMA_STREAM_2, LL_DMA_MEMORY_INCREMENT);

  LL_DMA_SetPeriphSize(DMA2, LL_DMA_STREAM_2, LL_DMA_PDATAALIGN_HALFWORD);

  LL_DMA_SetMemorySize(DMA2, LL_DMA_STREAM_2, LL_DMA_MDATAALIGN_HALFWORD);

  LL_DMA_DisableFifoMode(DMA2, LL_DMA_STREAM_2);

  /* USER CODE BEGIN TIM1_Init 1 */

  /* USER CODE END TIM1_Init 1 */
  TIM_InitStruct.Prescaler = 0;
  TIM_InitStruct.CounterMode = LL_TIM_COUNTERMODE_UP;
  TIM_InitStruct.Autoreload = 1904;
  TIM_InitStruct.ClockDivision = LL_TIM_CLOCKDIVISION_DIV1;
  TIM_InitStruct.RepetitionCounter = 0;
  LL_TIM_Init(TIM1, &TIM_InitStruct);
  LL_TIM_DisableARRPreload(TIM1);
  LL_TIM_OC_EnablePreload(TIM1, LL_TIM_CHANNEL_CH1);
  TIM_OC_InitStruct.OCMode = LL_TIM_OCMODE_PWM1;
  TIM_OC_InitStruct.OCState = LL_TIM_OCSTATE_DISABLE;
  TIM_OC_InitStruct.OCNState = LL_TIM_OCSTATE_DISABLE;
  TIM_OC_InitStruct.CompareValue = 0;
  TIM_OC_InitStruct.OCPolarity = LL_TIM_OCPOLARITY_HIGH;
  TIM_OC_InitStruct.OCNPolarity = LL_TIM_OCPOLARITY_HIGH;
  TIM_OC_InitStruct.OCIdleState = LL_TIM_OCIDLESTATE_LOW;
  TIM_OC_InitStruct.OCNIdleState = LL_TIM_OCIDLESTATE_LOW;
  LL_TIM_OC_Init(TIM1, LL_TIM_CHANNEL_CH1, &TIM_OC_InitStruct);
  LL_TIM_OC_DisableFast(TIM1, LL_TIM_CHANNEL_CH1);
  LL_TIM_OC_EnablePreload(TIM1, LL_TIM_CHANNEL_CH2);
  LL_TIM_OC_Init(TIM1, LL_TIM_CHANNEL_CH2, &TIM_OC_InitStruct);
  LL_TIM_OC_DisableFast(TIM1, LL_TIM_CHANNEL_CH2);
  LL_TIM_SetTriggerOutput(TIM1, LL_TIM_TRGO_RESET);
  LL_TIM_DisableMasterSlaveMode(TIM1);
  TIM_BDTRInitStruct.OSSRState = LL_TIM_OSSR_DISABLE;
  TIM_BDTRInitStruct.OSSIState = LL_TIM_OSSI_DISABLE;
  TIM_BDTRInitStruct.LockLevel = LL_TIM_LOCKLEVEL_OFF;
  TIM_BDTRInitStruct.DeadTime = 0;
  TIM_BDTRInitStruct.BreakState = LL_TIM_BREAK_DISABLE;
  TIM_BDTRInitStruct.BreakPolarity = LL_TIM_BREAK_POLARITY_HIGH;
  TIM_BDTRInitStruct.AutomaticOutput = LL_TIM_AUTOMATICOUTPUT_DISABLE;
  LL_TIM_BDTR_Init(TIM1, &TIM_BDTRInitStruct);
  /* USER CODE BEGIN TIM1_Init 2 */

  /* USER CODE END TIM1_Init 2 */
  LL_AHB1_GRP1_EnableClock(LL_AHB1_GRP1_PERIPH_GPIOB);
  LL_AHB1_GRP1_EnableClock(LL_AHB1_GRP1_PERIPH_GPIOA);
  /**TIM1 GPIO Configuration
  PB0   ------> TIM1_CH2N
  PB13   ------> TIM1_CH1N
  PA8   ------> TIM1_CH1
  PA9   ------> TIM1_CH2
  */
  GPIO_InitStruct.Pin = LL_GPIO_PIN_0|LL_GPIO_PIN_13;
  GPIO_InitStruct.Mode = LL_GPIO_MODE_ALTERNATE;
  GPIO_InitStruct.Speed = LL_GPIO_SPEED_FREQ_HIGH;
  GPIO_InitStruct.OutputType = LL_GPIO_OUTPUT_PUSHPULL;
  GPIO_InitStruct.Pull = LL_GPIO_PULL_NO;
  GPIO_InitStruct.Alternate = LL_GPIO_AF_1;
  LL_GPIO_Init(GPIOB, &GPIO_InitStruct);

  GPIO_InitStruct.Pin = LL_GPIO_PIN_8|LL_GPIO_PIN_9;
  GPIO_InitStruct.Mode = LL_GPIO_MODE_ALTERNATE;
  GPIO_InitStruct.Speed = LL_GPIO_SPEED_FREQ_HIGH;
  GPIO_InitStruct.OutputType = LL_GPIO_OUTPUT_OPENDRAIN;
  GPIO_InitStruct.Pull = LL_GPIO_PULL_NO;
  GPIO_InitStruct.Alternate = LL_GPIO_AF_1;
  LL_GPIO_Init(GPIOA, &GPIO_InitStruct);

}

/**
  * @brief TIM3 Initialization Function
  * @param None
  * @retval None
  */
static void MX_TIM3_Init(void)
{

  /* USER CODE BEGIN TIM3_Init 0 */

  /* USER CODE END TIM3_Init 0 */

  TIM_ClockConfigTypeDef sClockSourceConfig = {0};
  TIM_MasterConfigTypeDef sMasterConfig = {0};

  /* USER CODE BEGIN TIM3_Init 1 */

  /* USER CODE END TIM3_Init 1 */
  htim3.Instance = TIM3;
  htim3.Init.Prescaler = 210-1 ;
  htim3.Init.CounterMode = TIM_COUNTERMODE_UP;
  htim3.Init.Period = 65535;
  htim3.Init.ClockDivision = TIM_CLOCKDIVISION_DIV1;
  htim3.Init.AutoReloadPreload = TIM_AUTORELOAD_PRELOAD_DISABLE;
  if (HAL_TIM_Base_Init(&htim3) != HAL_OK)
  {
    Error_Handler();
  }
  sClockSourceConfig.ClockSource = TIM_CLOCKSOURCE_INTERNAL;
  if (HAL_TIM_ConfigClockSource(&htim3, &sClockSourceConfig) != HAL_OK)
  {
    Error_Handler();
  }
  sMasterConfig.MasterOutputTrigger = TIM_TRGO_RESET;
  sMasterConfig.MasterSlaveMode = TIM_MASTERSLAVEMODE_DISABLE;
  if (HAL_TIMEx_MasterConfigSynchronization(&htim3, &sMasterConfig) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN TIM3_Init 2 */

  /* USER CODE END TIM3_Init 2 */

}

/**
  * @brief USART1 Initialization Function
  * @param None
  * @retval None
  */
static void MX_USART1_UART_Init(void)
{

  /* USER CODE BEGIN USART1_Init 0 */

  /* USER CODE END USART1_Init 0 */

  /* USER CODE BEGIN USART1_Init 1 */

  /* USER CODE END USART1_Init 1 */
  huart1.Instance = USART1;
  huart1.Init.BaudRate = 115200;
  huart1.Init.WordLength = UART_WORDLENGTH_8B;
  huart1.Init.StopBits = UART_STOPBITS_1;
  huart1.Init.Parity = UART_PARITY_NONE;
  huart1.Init.Mode = UART_MODE_TX_RX;
  huart1.Init.HwFlowCtl = UART_HWCONTROL_NONE;
  huart1.Init.OverSampling = UART_OVERSAMPLING_16;
  if (HAL_UART_Init(&huart1) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN USART1_Init 2 */

  /* USER CODE END USART1_Init 2 */

}

/**
  * Enable DMA controller clock
  */
static void MX_DMA_Init(void)
{

  /* DMA controller clock enable */
  __HAL_RCC_DMA1_CLK_ENABLE();
  __HAL_RCC_DMA2_CLK_ENABLE();

  /* DMA interrupt init */
  /* DMA1_Stream4_IRQn interrupt configuration */
  HAL_NVIC_SetPriority(DMA1_Stream4_IRQn, 0, 0);
  HAL_NVIC_EnableIRQ(DMA1_Stream4_IRQn);
  /* DMA2_Stream1_IRQn interrupt configuration */
  NVIC_SetPriority(DMA2_Stream1_IRQn, NVIC_EncodePriority(NVIC_GetPriorityGrouping(),0, 0));
  NVIC_EnableIRQ(DMA2_Stream1_IRQn);
  /* DMA2_Stream2_IRQn interrupt configuration */
  NVIC_SetPriority(DMA2_Stream2_IRQn, NVIC_EncodePriority(NVIC_GetPriorityGrouping(),0, 0));
  NVIC_EnableIRQ(DMA2_Stream2_IRQn);

}

/**
  * @brief GPIO Initialization Function
  * @param None
  * @retval None
  */
static void MX_GPIO_Init(void)
{
  GPIO_InitTypeDef GPIO_InitStruct = {0};

  /* GPIO Ports Clock Enable */
  __HAL_RCC_GPIOC_CLK_ENABLE();
  __HAL_RCC_GPIOH_CLK_ENABLE();
  __HAL_RCC_GPIOA_CLK_ENABLE();
  __HAL_RCC_GPIOB_CLK_ENABLE();

  /*Configure GPIO pin Output Level */
  HAL_GPIO_WritePin(LED_GPIO_Port, LED_Pin, GPIO_PIN_RESET);

  /*Configure GPIO pin Output Level */
  HAL_GPIO_WritePin(GPIOA, LCD1_CMD_Pin|DUMMY_OUT_Pin, GPIO_PIN_RESET);

  /*Configure GPIO pin Output Level */
  HAL_GPIO_WritePin(MUTE_OUT_GPIO_Port, MUTE_OUT_Pin, GPIO_PIN_SET);

  /*Configure GPIO pin Output Level */
  HAL_GPIO_WritePin(GPIOB, LCD_CMD_Pin|LCD_RESET_Pin|LCD_BL_Pin, GPIO_PIN_RESET);

  /*Configure GPIO pin : LED_Pin */
  GPIO_InitStruct.Pin = LED_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_OD;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
  HAL_GPIO_Init(LED_GPIO_Port, &GPIO_InitStruct);

  /*Configure GPIO pin : KEY_Pin */
  GPIO_InitStruct.Pin = KEY_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_INPUT;
  GPIO_InitStruct.Pull = GPIO_PULLUP;
  HAL_GPIO_Init(KEY_GPIO_Port, &GPIO_InitStruct);

  /*Configure GPIO pin : LCD1_CMD_Pin */
  GPIO_InitStruct.Pin = LCD1_CMD_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_HIGH;
  HAL_GPIO_Init(LCD1_CMD_GPIO_Port, &GPIO_InitStruct);

  /*Configure GPIO pins : MUTE_OUT_Pin LCD_RESET_Pin LCD_BL_Pin */
  GPIO_InitStruct.Pin = MUTE_OUT_Pin|LCD_RESET_Pin|LCD_BL_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
  HAL_GPIO_Init(GPIOB, &GPIO_InitStruct);

  /*Configure GPIO pin : DUMMY_OUT_Pin */
  GPIO_InitStruct.Pin = DUMMY_OUT_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
  HAL_GPIO_Init(DUMMY_OUT_GPIO_Port, &GPIO_InitStruct);

  /*Configure GPIO pin : LCD_CMD_Pin */
  GPIO_InitStruct.Pin = LCD_CMD_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_HIGH;
  HAL_GPIO_Init(LCD_CMD_GPIO_Port, &GPIO_InitStruct);

}

/* USER CODE BEGIN 4 */

/* USER CODE END 4 */

/**
  * @brief  This function is executed in case of error occurrence.
  * @retval None
  */
void Error_Handler(void)
{
  /* USER CODE BEGIN Error_Handler_Debug */
  /* User can add his own implementation to report the HAL error return state */
  __disable_irq();
  while (1)
  {
  }
  /* USER CODE END Error_Handler_Debug */
}

#ifdef  USE_FULL_ASSERT
/**
  * @brief  Reports the name of the source file and the source line number
  *         where the assert_param error has occurred.
  * @param  file: pointer to the source file name
  * @param  line: assert_param error line source number
  * @retval None
  */
void assert_failed(uint8_t *file, uint32_t line)
{
  /* USER CODE BEGIN 6 */
  /* User can add his own implementation to report the file name and line number,
     ex: printf("Wrong parameters value: file %s on line %d\r\n", file, line) */
  /* USER CODE END 6 */
}
#endif /* USE_FULL_ASSERT */

/************************ (C) COPYRIGHT STMicroelectronics *****END OF FILE****/
