#include "ift.h"

/* it returns a set with external border pixels */
iftSet *MyBackgroundBorder(iftImage *bin)
{
    iftSet *BG_Set = NULL;
    iftAdjRel *A1 = iftCircular(1.0);
    int binSize = bin->n;

    for (int p = 0; p < binSize; p++)
    {
        if (bin->val[p] != 0)
        {
            for (int adj_idx = 0; adj_idx <= A1->n; adj_idx++)
            {
                // printf("ADJ\n");
                iftVoxel v_center = iftGetVoxelCoord(bin, p);
                iftVoxel adj_voxel = iftGetAdjacentVoxel(A1, v_center, adj_idx);
                // printf("got %d\n", adj_voxel.t);
                if (iftValidVoxel(bin, adj_voxel))
                {
                    int voxel_idx = iftGetVoxelIndex(bin, adj_voxel);
                    if (bin->val[voxel_idx] == 0)
                    {
                        // printf("Inserting\n");
                        iftInsertSet(&BG_Set, voxel_idx);
                    }
                }
            }
        }
    }

    return BG_Set;
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
        sprintf(filename, "%s/%s_select.png", "background", basename);
        iftWriteImageByExt(aux2, filename);

        iftSet *S_bg = NULL;
        S_bg = MyBackgroundBorder(aux2);
        aux2 = iftGrayImageToColorImage(aux2, iftGrayColorTable(256));
        for (iftSet *v_bg = S_bg; v_bg != NULL; v_bg = v_bg->next)
        {
            if (aux2->Cr == NULL) {
                printf("DAMN\n");
            }
            aux2->Cr[v_bg->elem] = 256;
        }
        sprintf(filename, "%s/%s_bg.png", "background", basename);
        iftWriteImageByExt(aux2, filename);

        iftDestroyImage(&aux1);
        iftDestroyImage(&aux2);
        iftDestroyImage(&orig);
        iftDestroyImage(&norm);
        iftFree(basename);

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
