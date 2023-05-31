
//-> in verschiedenen solvern vorgesehen; 

set_unknowns()
for(int i=0; i<unknownsDeltaU->MyLength(); i+=(3+numMultiphysDoFs)){
        for(int j=0; j<3; ++j){
          unknownsDeltaUPtr[i+j] = deltaUPtr[i*3/(3+numMultiphysDoFs)+j];
        }
        unknownsDeltaUPtr[i+3] = fluidPressureDeltaUPtr[i/(3+numMultiphysDoFs)];
      }

      for(int i=0 ; i<fluidPressureY->MyLength() ; ++i){
        fluidPressureYPtr[i] = fluidPressureUPtr[i] + fluidPressureDeltaUPtr[i];
        fluidPressureVPtr[i] = fluidPressureDeltaUPtr[i]/timeIncrement;
      }

      for(int i=0 ; i<unknownsY->MyLength() ; i+=(3+numMultiphysDoFs)){
        for(int j = 0; j<3; ++j){
          unknownsYPtr[i+j] = yPtr[i/(3+numMultiphysDoFs)*3 + j];
          unknownsVPtr[i+j] = vPtr[i/(3+numMultiphysDoFs)*3 + j];
        }
        unknownsYPtr[i+3] = fluidPressureYPtr[i/(3+numMultiphysDoFs)];
        unknownsVPtr[i+3] = fluidPressureVPtr[i/(3+numMultiphysDoFs)];
      }