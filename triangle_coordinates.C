
const int NLayer = 7;
const int NStaves[NLayer] = { 12, 16, 20, 24, 30, 42, 48 };
float MidPointRad[7] = { 23.49, 31.586, 39.341, 197.598, 246.944, 345.348, 394.883 };
const float StartAngle[7] = { 16.997 / 360 * (TMath::Pi() * 2.), 17.504 / 360 * (TMath::Pi() * 2.), 17.337 / 360 * (TMath::Pi() * 2.), 8.75 / 360 * (TMath::Pi() * 2.), 7 / 360 * (TMath::Pi() * 2.), 5.27 / 360 * (TMath::Pi() * 2.), 4.61 / 360 * (TMath::Pi() * 2.) };

void c(int layer, int stave, double* px, double* py)
{
  float stepAngle = TMath::Pi() * 2 / NStaves[layer];             // the angle between to stave
  float midAngle = StartAngle[layer] + (stave * stepAngle);       // mid point angle
  float staveRotateAngle = TMath::Pi() / 2 - (stave * stepAngle); // how many angle this stave rotate(compare with first stave)
  px[1] = MidPointRad[layer] * TMath::Cos(midAngle);              // there are 4 point to decide this TH2Poly bin
                                                                  // 0:left point in this stave;
                                                                  // 1:mid point in this stave;
                                                                  // 2:right point in this stave;
                                                                  // 3:higher point int this stave;
  py[1] = MidPointRad[layer] * TMath::Sin(midAngle);              // 4 point calculated accord the blueprint
                                                                  // roughly calculate
  if (layer < 3) {
    px[0] = 7.7 * TMath::Cos(staveRotateAngle) + px[1];
    py[0] = -7.7 * TMath::Sin(staveRotateAngle) + py[1];
    px[2] = -7.7 * TMath::Cos(staveRotateAngle) + px[1];
    py[2] = 7.7 * TMath::Sin(staveRotateAngle) + py[1];
    px[3] = 5.623 * TMath::Sin(staveRotateAngle) + px[1];
    py[3] = 5.623 * TMath::Cos(staveRotateAngle) + py[1];
  } else {
    px[0] = 21 * TMath::Cos(staveRotateAngle) + px[1];
    py[0] = -21 * TMath::Sin(staveRotateAngle) + py[1];
    px[2] = -21 * TMath::Cos(staveRotateAngle) + px[1];
    py[2] = 21 * TMath::Sin(staveRotateAngle) + py[1];
    px[3] = 40 * TMath::Sin(staveRotateAngle) + px[1];
    py[3] = 40 * TMath::Cos(staveRotateAngle) + py[1];
  }
}


void triangle_coordinates(){

  int counter = 0;
    for (int ilayer = 0; ilayer < 7; ilayer++) {
      for (int istave = 0; istave < NStaves[ilayer]; istave++) {
	double *px = new double[4];
	double *py = new double[4];
	c(ilayer, istave, px, py);
	if (ilayer < 3) {
          for (int icoo = 0; icoo < 4; icoo++) {
            px[icoo] *= 3.;
            py[icoo] *= 3.;
          }
        }
	// point [1] is the middle point of the stave: not necessary
	cout<<counter<<" "<<px[0]<<" "<<px[2]<<" "<<px[3]<<" "<<py[0]<<" "<<py[2]<<" "<<py[3]<<endl;
	counter++;
      }
    }

    
	

	
}

