#include "rms.h"

/* CIRCULAR BUFFER FUNCTIONS IMPLEMENTATION */

int circularBufferInit(CircularBuffer *cb, int size)
{
  cb->buff = (float *)malloc(sizeof(float) * size);
  if (sizeof(cb->buff) == 0)
    return -1;
  cb->size = size;
  cb->start = 0;
  cb->end = 0;
  cb->count = 0;
  return 0;
}

int circularBufferIsEmpty(const CircularBuffer *cb)
{
  return cb->size == 0;
}

int circularBufferIsFull(const CircularBuffer *cb)
{
  return cb->count == cb->size;
}

void add2CircularBuffer(CircularBuffer *cb, float data)
{
  cb->buff[cb->end] = data;
  cb->end = (cb->end + 1) % cb->size;
  if (cb->count == cb->size)
    cb->start = (cb->start + 1) % cb->size;
  else
    cb->count++;
}

float getElement(CircularBuffer *cb, int n)
{
  if (n >= cb->count || n < 0)
    return 0xFFFF;
  int index = (cb->start + n) % cb->size;
  return cb->buff[index];
}

void circularBufferFree(CircularBuffer *cb)
{
  free(cb->buff);
}

void plotBufferFromStartToEnd(CircularBuffer *cb)
{
  int i;
  int count = cb->count;
  if (count == 0)
  {
    printf("Buffer is empty.\n");
    return;
  }
  printf("Buffer content:\n");
  for (i = cb->start; count > 0; i = (i + 1) % cb->size)
  {
    printf("%d ", cb->buff[i]);
    count--;
  }
  printf("\n");
}

/* RMS CALCULATOR FUNCTIONS IMPLEMENTATION */

int RMSCalculatorInit(RMSCalculator *rms, int size)
{
  if (circularBufferInit(&rms->square_buffer, size))
    return -1;
  rms->size = size;
  rms->sum = 0;
}

void add2RMSCalculator(RMSCalculator *rms, float data)
{
  float value = data * data;
  float sub_value = 0;

  if (circularBufferIsFull(&rms->square_buffer))
    sub_value = getElement(&rms->square_buffer, 0);

  add2CircularBuffer(&rms->square_buffer, value);

  rms->sum = rms->sum + value - sub_value;
}

void calculateRMSvalue(RMSCalculator *rms, float *rms_value)
{
  float value = rms->sum / rms->size;
  *rms_value = sqrt(value);
}

void RMSCalculatorFree(RMSCalculator *rms)
{
  circularBufferFree(&rms->square_buffer);
}