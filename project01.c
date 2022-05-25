#include "ift.h"

/*
   Project 01, course MO445, Prof. Alexandre Falcao.

   This code crops a region of interest containing a fingerprint with
   as minimum as possible surrounding noise.

   The code contains three functions iftAsfCOBin(), iftCloseBasins(),
   iftErodeBin() that are implemented by the Image Foresting
   Transform. Your task is to substitute them by the corresponding
   functions implemented by you. In order to do that, you should fill
   the code of the functions below. You should use iftAddFrame for
   padding zeroes and iftRemFrame to return the original image size
   whenever needed. Object pixels in these functions are pixels with
   value different from zero (e.g., usually 1 or 255) and background
   pixels are those with value equal to zero.  Object pixels are
   internal border pixels when they have a four-neighbor outside the
   object and background pixels are external border pixels when they
   have a four-neighbor inside an object.
 */

/* it returns pixels at the border of the image */
iftSet *MyImageBorder(iftImage *bin)
{
   iftSet *imgBorder = NULL;
   return imgBorder;
}


/* it returns a set with external border pixels */
iftSet *MyBackgroundBorder(iftImage *bin)
{
   iftSet *objBorder = NULL;
   return objBorder;
}

/* it returns a set with internal border pixels */
iftSet *MyObjectBorder(iftImage *bin)
{
   iftSet *S = NULL;
   iftAdjRel *A1 = iftCircular(1.0);
   iftVoxel p_voxel, q_voxel;
   int binSize = bin->n;

   for (int p = 0; p < binSize; p++)
   {
      if (bin->val[p] == 0)
         continue;
      for (int adj_idx = 0; adj_idx <= A1->n; adj_idx++)
      {
         p_voxel = iftGetVoxelCoord(bin, p);
         q_voxel = iftGetAdjacentVoxel(A1, p_voxel, adj_idx);
         if (iftValidVoxel(bin, q_voxel))
         {
            int q_idx = iftGetVoxelIndex(bin, q_voxel);
            int p_idx = iftGetVoxelIndex(bin, p_voxel);
            // Adjecent pixel is background
            if (bin->val[q_idx] == 0)
            {
               iftInsertSet(&S, p_idx);
            }
         }
      }
   }

   iftDestroyAdjRel(&A1);

   return S;
}

iftImage *MyDilateBin(iftImage *bin, iftSet **S, float radius)
{
   iftImage *D = iftCopyImage(bin);
   int sqr_radius = pow(radius,2);
   iftAdjRel *Asqrt = iftCircular(sqrt(2));
   iftImage *C = NULL;
   iftImage *R = NULL;
   iftGQueue *Q = NULL;
   C = iftCreateImage(bin->xsize, bin->ysize, bin->zsize);
   R = iftCreateImage(bin->xsize, bin->ysize, bin->zsize);
   // Queue parameters
   // Priority queue
   Q = iftCreateGQueue(bin->n + 1, bin->n, C->val);

   // Get External border
   *S = MyObjectBorder(bin);
   
   // Inint costs and root maps
   for (int p = 0; p < bin->n; p++)
   {
      C->val[p] = INT_MAX;
      R->val[p] = p;
   }
   
   // Set cost = 0 &
   //     root as itself
   // for internal border pixels
   // All min costs at queue start with C=0
   int p = 0;
   while (*S != NULL)
   {
      p = iftRemoveSet(S);
      C->val[p] = 0;
      R->val[p] = p;
      iftInsertGQueue(&Q, p);
   }
   
   // Dilate
   while (!iftEmptyGQueue(Q))
   {
      int p = iftRemoveGQueue(Q);
      if (C->val[p] <= sqr_radius)
      {
         D->val[p] = bin->val[R->val[p]];

         for (int adj_idx = 0; adj_idx <= Asqrt->n; adj_idx++)
         {
            iftVoxel p_voxel = iftGetVoxelCoord(D, p);
            iftVoxel q_voxel = iftGetAdjacentVoxel(Asqrt, p_voxel, adj_idx);
            
            if (!iftValidVoxel(D, q_voxel))
               continue;
            
            int q = iftGetVoxelIndex(D, q_voxel);
            if (C->val[q] > C->val[p] && bin->val[q] == 0)
            {
               iftVoxel rp = iftGetVoxelCoord(D, R->val[p]);
               int tmp = pow(p_voxel.x - rp.x, 2) + pow(p_voxel.y - rp.y, 2);
               if (tmp < C->val[q])
               {
                  if (Q->L.elem[0].color == IFT_GRAY)
                     iftRemoveGQueue(Q);
                  C->val[q] = tmp;
                  R->val[q] = R->val[p];
                  iftInsertGQueue(&Q, q);
               }
            }
         }
      }
      else
      {
         iftInsertSet(S, p);
      }
   }

   // sprintf(filename, "dilation/%s_dilated.png", out_dir, basename);
   iftWriteImageByExt(D, "dilation/dilated.png");
   iftWriteImageByExt(bin, "dilation/non_dilated.png");

   free(C);
   free(R);
   iftDestroyIntQueue(&Q);
   exit(-1);

   return D;
}

/* it erodes objects */
iftImage *MyErodeBin(iftImage *bin, iftSet **S, float radius)
{
   iftImage *erodedBin = iftCopyImage(bin);
   return erodedBin;
}

/* it executes dilation followed by erosion */
iftImage *MyCloseBin(iftImage *bin, float radius)
{
   iftImage *dilatedBin, *closedBin;
   iftSet *S = NULL;

   dilatedBin = MyDilateBin(bin, &S, radius);
   closedBin = MyErodeBin(dilatedBin, &S, radius);

   free(dilatedBin);
   iftDestroySet(&S);

   return closedBin;
}

/* it executes erosion followed by dilation */
iftImage *MyOpenBin(iftImage *bin, float radius)
{
   iftImage *erodedBin, *openBin;
   iftSet *S = NULL;

   erodedBin = MyErodeBin(bin, &S, radius);
   openBin = MyDilateBin(erodedBin, &S, radius);

   free(erodedBin);
   iftDestroySet(&S);

   return openBin;
}

/* it executes closing followed by opening */
iftImage *MyAsfCOBin(iftImage *bin, float radius)
{
   iftImage *closedBin, *asfCOBin;

   closedBin = MyCloseBin(bin, radius);
   asfCOBin = MyOpenBin(closedBin, radius);

   free(closedBin);

   return asfCOBin;
}

/* it closes holes in objects */
iftImage *MyCloseBasins(iftImage *bin)
{
   iftImage *closedBasins = iftCopyImage(bin);
   return closedBasins;
}

int main(int argc, char *argv[])
{
   timer *tstart = NULL;
   char filename[200];

   /*--------------------------------------------------------*/

   void *trash = malloc(1);
   struct mallinfo info;
   int MemDinInicial, MemDinFinal;
   free(trash);
   info = mallinfo();
   MemDinInicial = info.uordblks;

   /*--------------------------------------------------------*/

   if (argc != 3)
   {
      printf("project01 <P1> <P2>\n");
      printf("P1: folder with original images\n");
      printf("P2: folder with cropped images\n");
      exit(0);
   }

   tstart = iftTic();

   iftFileSet *fs = iftLoadFileSetFromDirBySuffix(argv[1], ".png", 1);
   int nimages = fs->n;
   char *out_dir = argv[2];
   iftMakeDir(out_dir);
   iftAdjRel *A = iftCircular(3.5), *B = iftCircular(1.5);

   for (int i = 0; i < nimages; i++)
   {
      char *basename = iftFilename(fs->files[i]->path, ".png");
      iftImage *orig = iftReadImageByExt(fs->files[i]->path);
      /* normalize  image */
      iftImage *norm = iftNormalize(orig, 0, 255);
      /* binarize image */
      iftImage *aux1 = iftBelowAdaptiveThreshold(norm, NULL, A, 0.98, 2, 255);
      /* remove noise components from the background */
      iftImage *aux2 = iftSelectCompAboveArea(aux1, B, 100);
      iftDestroyImage(&aux1);
      // sprintf(filename, "%s/%s_select.png", out_dir, basename);
      // iftWriteImageByExt(aux2, filename);

      /* apply morphological filtering to make the fingerprint the
         largest component: this operation must add frame and remove it
         afterwards. */
      aux1 = MyAsfCOBin(aux2, 15.0); // iftAsfCOBin(aux2,15.0);
      iftDestroyImage(&aux2);
      // sprintf(filename, "%s/%s_asfCOBin.png", out_dir, basename);
      // iftWriteImageByExt(aux1, filename);
      /* close holes inside the components to allow subsequent erosion
         from the external borders only */
      aux2 = MyCloseBasins(aux1); // iftCloseBasins(aux1,NULL,NULL);
      iftDestroyImage(&aux1);
      // sprintf(filename, "%s/%s_CloseBasins.png", out_dir, basename);
      // iftWriteImageByExt(aux2, filename);
      /* erode components and select the largest one to estimate its
         center as close as possible to the center of the fingerprint */
      iftSet *S = NULL;
      aux1 = MyErodeBin(aux2, &S, 30.0); // iftErodeBin(aux2,&S,30.0);

      // sprintf(filename, "%s/%s_ErodeBin.png", out_dir, basename);
      // iftWriteImageByExt(aux1, filename);
      iftDestroySet(&S);
      iftDestroyImage(&aux2);
      aux2 = iftSelectLargestComp(aux1, B);

      /* crop the normalized image by the minimum bounding box of the
         resulting mask (largest component) */

      iftDestroyImage(&aux1);
      iftVoxel pos;
      iftBoundingBox bb = iftMinBoundingBox(aux2, &pos);
      aux1 = iftExtractROI(norm, bb);

      sprintf(filename, "%s/%s.png", out_dir, basename);
      iftWriteImageByExt(aux1, filename);
      iftDestroyImage(&aux1);
      iftDestroyImage(&aux2);
      iftDestroyImage(&orig);
      iftDestroyImage(&norm);
      iftFree(basename);
      printf("Done %d images.\n", i + 1);
      break;
   }

   iftDestroyFileSet(&fs);
   iftDestroyAdjRel(&A);
   iftDestroyAdjRel(&B);

   puts("\nDone...");
   puts(iftFormattedTime(iftCompTime(tstart, iftToc())));

   /* ---------------------------------------------------------- */

   info = mallinfo();
   MemDinFinal = info.uordblks;
   if (MemDinInicial != MemDinFinal)
      printf("\n\nDinamic memory was not completely deallocated (%d, %d)\n",
             MemDinInicial, MemDinFinal);

   return 0;
}
