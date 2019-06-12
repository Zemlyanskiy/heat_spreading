#pragma once
#include "IOStream.h"
#include "MainVaribles.h"

double CalcFunc(double Xcoord, double Ycoord, double Zcoord)
{
  if (Xcoord >= -0.5 && Xcoord <= 0.5 &&
      Ycoord >= -0.5 && Ycoord <= 0.5 &&
      Zcoord >= -0.5 && Zcoord <= 0.5)
  {
    return cos(Xcoord * 3.141592) + cos(Ycoord * 3.141592) + cos(Zcoord * 3.141592);
  }
  return 0;
}

static inline void StartString(struct input_data temp)
{
    posX = temp.XMin;
    posY = temp.YMin;
    posZ = temp.ZMin;
    for (int i = 0; i < temp.array.Lx; i++)
    {
        posX += Xstep;
        posY = temp.YMin;
        for (int j = 0; j < temp.array.Ly; j++)
        {
            posY += Ystep;
            posZ = temp.ZMin;
            for (int k = 0; k < temp.array.Lz; k++)
            {
                posZ += Zstep;
                SetValue(&points, i, j, k, CalcFunc(posX, posY, posZ));
            }
        }
    }
    //First output result
    Output(file, points.values, temp.array.Lx * temp.array.Ly * temp.array.Lz);
}
