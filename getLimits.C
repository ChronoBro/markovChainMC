void getLimits(){

  ostringstream inputFileName;
  inputFileName << "outputs/doubletGroundState/reduceTest/probDensity0.txt";
  //inputFileName << "outputs/tripletGroundState/reduceTest/probDensity0.txt";

  ifstream inputFile(inputFileName.str().c_str());

  double x[100];
  double y[100];

 
  int i=0;
  while(!inputFile.eof()){
    inputFile >> x[i] >> y[i];
    i++;

  }

  TGraph * gr = new TGraph(100,x,y);

  gr->Draw("A*");

}
