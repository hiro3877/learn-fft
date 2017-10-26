//
//  bitmap structs
//
// (c)Copyright Spacesoft corp., 2017 All rights reserved.
//                               Hiro KITAYAMA
#ifndef __BITMAPSTRUCT__
#define __BITMAPSTRUCT__

#pragma pack(push, 1)

typedef struct
{
    unsigned short  bfType;
    unsigned int    bfSize;
    unsigned short  bfReserved1;
    unsigned short  bfReserved2;
    unsigned int    bfOffBits;
}
bmpFileHdr, *pBmpFileHdr;

typedef struct
{
    unsigned int    biSize;
    int             biWidth;
    int             biHeight;
    unsigned short  biPlanes;
    unsigned short  biBitCount;
    unsigned int    biCompression;
    unsigned int    biSizeImage;
    int             biXPelsPerMeter;
    int             biYPelsPerMeter;
    unsigned int    biClrUsed;
    unsigned int    biClrImportant;
}
bmpInfoHdr, *pBmpInfoHdr;

#pragma pack(pop)

#define SAFE_DELETE(p)    if(p!=0) { delete(p); p=0; }
#define SAFE_FREE(p)      if(p!=0) { free(p);   p=0; }

#endif // __BITMAPSTRUCT__
