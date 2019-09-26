void filter(){

  ostringstream name;
  name << "decayEnergy_1keVBins_filter.txt";
  ifstream infile(name.str().c_str());
  ostringstream name2;
  name2 << "decayEnergy_1keVBins_filter2.txt";
  ofstream ofile(name2.str().c_str());
  
  for(int i=0;!infile.eof();i++){

    double x,y;
    infile >> x;
    infile >> y;

    if(x<800)
      continue;

    if(x>1800)
      break;
    
    if(y>0)
      ofile << x << "\t" << y << "\n";
    
    
  }
  
}
