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
   iftImage *imgBorder;
   return imgBorder;
}

/* it returns a set with internal border pixels */
iftSet *MyObjectBorder(iftImage *bin)
{
   iftImage *objBorder;
   return objBorder;
}

/* it returns a set with external border pixels */
iftSet *MyBackgroundBorder(iftImage *bin)
{
   iftImage *bgBorder;
   return bgBorder;
}

/* it dilates objects */
iftImage *MyDilateBin(iftImage *bin, iftSet **S, float radius)
{
   iftImage *dilatedBin = iftCopyImage(bin);
   return dilatedBin;
}

/* it erodes objects */
iftImage *MyErodeBin(iftImage *bin, iftSet **S, float radius)
{
   iftImage *erodedBin;
   return erodedBin;
}

/* it executes dilation followed by erosion */
iftImage *MyCloseBin(iftImage *bin, float radius)
{
   iftImage *closedBin;
   return closedBin;
}

/* it executes erosion followed by dilation */
iftImage *MyOpenBin(iftImage *bin, float radius)
{
   iftImage *openBin;
   return openBin;
}

/* it executes closing followed by opening */
iftImage *MyAsfCOBin(iftImage *bin, float radius)
{
   iftImage *asfBin;
   return asfBin;
}

/* it closes holes in objects */
iftImage *MyCloseBasins(iftImage *bin)
{
   iftImage *closedBasins;
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
      sprintf(filename, "%s/%s_select.png", out_dir, basename);
      iftWriteImageByExt(aux2, filename);

      /* apply morphological filtering to make the fingerprint the
         largest component: this operation must add frame and remove it
         afterwards. */
      aux1 = MyAsfCOBin(aux2, 15.0); // iftAsfCOBin(aux2,15.0);
      iftDestroyImage(&aux2);
      /* close holes inside the components to allow subsequent erosion
         from the external borders only */
      aux2 = MyCloseBasins(aux1); // iftCloseBasins(aux1,NULL,NULL);
      iftDestroyImage(&aux1);
      /* erode components and select the largest one to estimate its
         center as close as possible to the center of the fingerprint */
      iftSet *S = NULL;
      aux1 = MyErodeBin(aux2, &S, 30.0); // iftErodeBin(aux2,&S,30.0);
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
