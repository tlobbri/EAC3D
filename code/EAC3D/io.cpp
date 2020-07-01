
#include "stdafx.h"


#include "EAC3D.h"
#include <string.h>
#include <iostream>


//using namespace std;
size_t fwriteStr80(char *str, FILE *fp){

  char cbuf[81];
  size_t len;
  int i;

  /* Most of this just to avoid garbage after the name. */
  len = strlen(str);
  strncpy(cbuf, str, len < 80 ? len : 80);

  for (i = len; i < 80; i++){cbuf[i] = '\0';}  /* pad with zeros */

  return fwrite(cbuf, sizeof(char), 80, fp);

}/* end fwriteStr80 */


void EAC3D::saveMesh(){
  int nbBloc = 1;
  int ibuf[3] ;
  float fbuf[3] ;
  char sbuf[81] ;

  char fileName[80] ;
  FILE *fmesh ;

    int id, jd, kd,tmp;
   
    sprintf(fileName, "mesh.dat");


    fmesh = fopen(fileName,"wb") ;

    /* write Ensight Geo Header */

    fwriteStr80((char *)"C Binary",fmesh);
    fwriteStr80((char *)"ENSIGHT GOLD GEOMETRY file format: CASE binary fileformat\n",fmesh);
    sprintf(sbuf,"EAC3D cartesian meshes" );
    fwriteStr80(sbuf,fmesh);
    fwriteStr80((char *)"node id assign",fmesh);
    fwriteStr80((char *)"element id assign",fmesh);

  for (int n=0; n<nbBloc; n++){
   
      /*##########################################################################*/
       /* Writing part metric */
      fwriteStr80((char *)"part",fmesh);
      /* Part id number*/
      ibuf[0] = n+1;
      fwrite(ibuf, sizeof(int), 1, fmesh);
      /* Comment : part name */
      sprintf(sbuf,"mesh_%d", n);
      fwriteStr80(sbuf,fmesh);
      /* Block type */
      fwriteStr80((char *)"block rectilinear", fmesh); //with_ghost

      /* Block size */
      ibuf[0] = this->nx;
      ibuf[1] = this->ny;//ngsize;//GNGsize;
      ibuf[2] = this->nz;//ogsize;//GOGsize;
      fwrite(ibuf, sizeof(int), 3, fmesh);

      
      // Block x
	  for (int i = 0; i < this->nx ; i++) {
		  fbuf[0] = i*dx + this->x0;
		  fwrite(fbuf, sizeof(float), 1, fmesh);
	  }

	  //Block y
	  for (int j = 0; j < this->ny ; j++) {
		  fbuf[0] = j*dy + this->y0;
		  fwrite(fbuf, sizeof(float), 1, fmesh);
	  }
	  //Block z
	  for (int k = 0; k < this->nz ; k++) {
		  fbuf[0] = k *dz + this->z0;
		  fwrite(fbuf, sizeof(float), 1, fmesh);
	  }

	  
  }
  fclose(fmesh);

}


void EAC3D::saveData(std::string s, int bnd_id) {
	std::cout << s << std::endl;

	FILE *fmesh;
	const char* fileName = s.c_str();
	fmesh = fopen(fileName, "wb");
	fwriteStr80((char *)"EAC3D : EnSight Gold scalar file format", fmesh);
	int nbBloc = 1;

	int ibuf[3];
	float fbuf[3];

	/* Writing part metric */
	for (int n = 0; n < nbBloc; n++) {

		/* part header */
		fwriteStr80((char *)"part", fmesh);

		*ibuf = n + 1;
		fwrite(ibuf, sizeof(int), 1, fmesh);
		fwriteStr80((char *)"block", fmesh);

		for (int k = 0; k < this->nz; k++) {
			for (int j = 0; j < this->ny; j++) {
				for (int i = 0; i < this->nx; i++) {

					if (bnd_id == STRESSID) {
						*fbuf = (float)(this->Stress[i][j][k].back());
					}
					else if (bnd_id == STRAINID) {
						*fbuf = (float)(this->Strain[i][j][k].back());
					}
					else {
						*fbuf = (float)(this->F[bnd_id][i][j][k]);
					}
					fwrite(fbuf, sizeof(float), 1, fmesh);
				}
			}
		}
	}
	fclose(fmesh);


}


