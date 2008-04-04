/* -*- mode: c++ -*- */

#include "RecoBTag/PerformanceMeasurements/interface/TtTagConsistencyEfficiency.h"

TtTagConsistencyEfficiency::_tableLine::_tableLine()
{
  discr = 0.0;
  for ( int i = 0; i < 5; i++ ) N[i] = 0;
  Nb = 0;
  Nbtag = 0;
  Nc = 0;
  Nctag = 0;
  Nl = 0;
  Nltag = 0;
  Nx = 0;
  Nxtag = 0;
}


TtTagConsistencyEfficiency::Fijk::Fijk()
{
  nEntries = 0;
  sum = 0;
}


void TtTagConsistencyEfficiency::cleanTable( tableContent & tab )
{
  for (int i = 0; i < 5; i++ )
    {
      tab . N[i] = 0;
    }
  tab . Nb = 0;
  tab . Nc = 0;
  tab . Nl = 0;
  tab . Nbtag = 0;
  tab . Nctag = 0;
  tab . Nltag = 0;
}

void TtTagConsistencyEfficiency::cleanTable( tableLine & tab )
{
  for (int i = 0; i < 5; i++ )
    {
      tab . N[i] = 0;
    }
  tab . Nb = 0;
  tab . Nc = 0;
  tab . Nl = 0;
  tab . Nx = 0;
  tab . Nbtag = 0;
  tab . Nctag = 0;
  tab . Nltag = 0;
  tab . Nxtag = 0;
}

TtTagConsistencyEfficiency::TtTagConsistencyEfficiency( void )
{
  inputFileRead = false;
}

int TtTagConsistencyEfficiency::test( void )
{
  for (int i = 0; i < Fijk_bg . nEntries; i++ )
    {
      cout << "F_" << Fijk_bg.i[i] << "_" << Fijk_bg.j[i] << "_" << Fijk_bg.k[i] << " = ";
      cout << Fijk_bg.value[Fijk_bg.i[i]][Fijk_bg.j[i]][Fijk_bg.k[i]] << endl;
    }

  return 0;
}

TtTagConsistencyEfficiency::BFormula TtTagConsistencyEfficiency::getFormula( int i, int j, int k, const char * suffix )
{
  ostringstream name, title, formula[5];
  name . str(""); // clear the buffer
  title . str("");
  for (int tagged = 0; tagged < 5; tagged++) formula[tagged] . str("");

  char buf[3];
  sprintf(buf, "%i%i%i", i,j,k);
  string Fijk = "F";
  Fijk = Fijk + buf + suffix;

  name << Fijk;
  title << Fijk;

  for (int tagged = 0; tagged < 5; tagged++)
    {

      formula[tagged] << Fijk << "*( ";
      
      bool firstEntry = true;
      
      for ( int ii=0; ii <= min(i,tagged); ii++ )
	{
	  for ( int jj=0; jj <= min(j,tagged-ii); jj++)
	    {
	      if ( (tagged - ii - jj) <= k )
		{
		  
		  int kk = tagged - ii - jj;
		  
		  int coeff = (int)( TMath::Binomial(i,ii) * TMath::Binomial(j,jj) * TMath::Binomial(k,kk) );
		  
		  if ( !firstEntry )
		    {
		      formula[tagged] << " + ";
		    }
		  
		  firstEntry = false;
		  
		  formula[tagged] << coeff;
		  
		  if (ii == 1)
		    {
		      formula[tagged] << " * epsb";
		    }
		  else if (ii>0)
		    {
		      formula[tagged] << " * epsb^" << ii;
		    }
		  
		  if (i-ii == 1)
		    {
		      formula[tagged] << " * (1-epsb)";
		    }
		  else if (i-ii>0)
		    {
		      formula[tagged] << " * (1-epsb)^" << (i-ii);
		    }
		  
		  if (jj == 1)
		    {
		      formula[tagged] << " * epsc";
		    }
		  else if (jj>0)
		    {
		      formula[tagged] << " * epsc^" << jj;
		    }
		  
		  if (j-jj == 1)
		    {
		      formula[tagged] << " * (1-epsc)";
		    }
		  else if (j-jj>0)
		    {
		      formula[tagged] << " * (1-epsc)^" << (j-jj);
		    }
		  
		  if (kk == 1)
		    {
		      formula[tagged] << " * epsl";
		    }
		  else if (kk>0)
		    {
		      formula[tagged] << " * epsl^" << kk;
		    }
		  
		  if (k-kk == 1)
		    {
		      formula[tagged] << " * (1-epsl)";
		    }
		  else if (k-kk>0)
		    {
		      formula[tagged] << " * (1-epsl)^" << (k-kk) ;
		    }
		  
		}
	    }
	}
      
      formula[tagged] << ")";
      
    }
  BFormula theFormula;
  theFormula . name = name . str();
  theFormula . title = title . str();
  for ( int ii = 0; ii < 5; ii++ ) theFormula . formula[ii] = formula[ii] . str();

  //for ( int l = 1; l < 4; l++ )
  //{
  //  cout << theFormula . name . c_str() << " #-----# " << theFormula . formula[l] . c_str() << endl;
  //}
  
  return theFormula;
  //return 1;
}

int TtTagConsistencyEfficiency::updateFijk( string mcListFileName, double weight, string dataType, int minValue )
{
  int result = 0;

  string buf;
  struct Fijk * _fijk = NULL;
  if ( dataType == "signal" ) _fijk = &Fijk_sig;
  else if ( dataType == "background" ) _fijk = &Fijk_bg;
  else
    {
      cout << "Unknown data type, exiting..." << endl;
      result = -1;
    }
  if ( _fijk )
    {
      int i,j,k;
      int value = 0;
      const char * _format = "F_%d_%d_%d = %d";

      // read the .mc file names from the list
      vector<string> mcFileName;
      char filename[1024];
      ifstream inFile( mcListFileName . c_str(), ios::in );
      if (!inFile)
	{
	  cout << " Unable to open list file: " << mcListFileName << endl;
	  return kFALSE;
	}
      else
	{
	  cout << "List file opened successfully: " << mcListFileName << endl;
	}
      while (inFile >> filename)
	{
	  string fullFileName = "";
	  fullFileName . append( filename );
	  mcFileName . push_back( fullFileName );
	}
      inFile.close();

      // reading files with Fijk
      vector<string>::const_iterator _filename;
      for ( _filename = mcFileName . begin(); _filename != mcFileName . end(); _filename++  )
	{
	  ifstream _mcFile( _filename -> c_str() );
	    
	  if ( _mcFile . is_open() )
	    {
	      
	      while ( getline( _mcFile, buf ) > 0 )
		{
		  sscanf( buf . c_str(), _format, &i, &j, &k, &value );
		  //TEMPORARY!!!
		  //if ( value >= minValue )
		  if ( value >= minValue && (i+j+k) >= 4)
		    {
		      int _index = -1;
		      for ( int _i = 0; (_i < _fijk -> nEntries) && _index < 0; _i++ )
			{
			  if ( (_fijk->i[_i] == i) &&(_fijk->j[_i] == j) && (_fijk->k[_i] == k) ) _index = _i;
			}
		      if ( _index < 0 )
			{
			  _fijk -> i . push_back( i );
			  _fijk -> j . push_back( j );
			  _fijk -> k . push_back( k );
			  _fijk -> value[i][j][k] = (int)( (double)value * weight );
			  _fijk -> nEntries++;
			}
		      else
			{
			  _fijk -> value[i][j][k] += (int)( (double)value * weight );
			}
		      _fijk -> sum += (int)( (double)value * weight );
		    }
		}
	      
	      cout << "Number of Fijk read: " << _fijk -> nEntries << endl;
	      
	      result = _fijk -> nEntries;
	    }
	}
    }
  
  return result;
}

int TtTagConsistencyEfficiency::updateFijk2( string mcListFileName, double weight, string dataType, int minValue )
{
  int result = 0;
  string buf;
  struct Fijk * _fijk = NULL;
  if ( dataType == "signal" ){
    _fijk = &Fijk_sig;
  }
  else if ( dataType == "background" ) _fijk = &Fijk_bg;
  else
    {
      cout << "Unknown data type, exiting..." << endl;
      result = -1;
    }
  if ( _fijk )
    {
      int i,j,k,x;
      int value = 0;
      const char * _format = "F_%d_%d_%d_%d = %d";

      // read the .mc file names from the list
      vector<string> mcFileName;
      char filename[1024];
      ifstream inFile( mcListFileName . c_str(), ios::in );
      if (!inFile)
	{
	  cout << " Unable to open list file: " << mcListFileName << endl;
	  return kFALSE;
	}
      else
	{
	  cout << "List file opened successfully: " << mcListFileName << endl;
	}
      while (inFile >> filename)
	{
	  string fullFileName = "";
	  fullFileName . append( filename );
	  mcFileName . push_back( fullFileName );
	}
      inFile.close();

      // reading files with Fijk
      vector<string>::const_iterator _filename;
      for ( _filename = mcFileName . begin(); _filename != mcFileName . end(); _filename++  )
	{

	  // quick fix to use .mc2 files from .mc file lists
	  string mc2FileName = *(_filename);
	  //mc2FileName . append( "2" ); // use file.mc2 instead of file.mc, it has unidentified jets

	  //ifstream _mcFile( _filename -> c_str() );
	  ifstream _mcFile( mc2FileName . c_str() );
	    
	  if ( _mcFile . is_open() )
	    {
	      while ( getline( _mcFile, buf ) > 0 )
		{
		  sscanf( buf . c_str(), _format, &x, &i, &j, &k, &value );
		  //TEMPORARY... later addition: ...+x must be ok
		  //if ( value >= minValue )
		  if ( value >= minValue && (i+j+k+x) >= 4)
		    {
		      int _index = -1;
		      for ( int _i = 0; (_i < _fijk -> nEntries) && _index < 0; _i++ )
			{
			  if ( (_fijk->i[_i] == i) &&(_fijk->j[_i] == j) && (_fijk->k[_i] == k+x) ) _index = _i;
			}
		      if ( _index < 0 )
			{
			  _fijk -> i . push_back( i );
			  _fijk -> j . push_back( j );
			  _fijk -> k . push_back( k+x );
			  _fijk -> value[i][j][k+x] = (int)( (double)value * weight );
			  _fijk -> nEntries++;
			}
		      else
			{
			  _fijk -> value[i][j][k+x] += (int)( (double)value * weight );
			}
		      _fijk -> sum += (int)( (double)value * weight );
		    }
		}
	      
	      cout << "Number of Fijk read: " << _fijk -> nEntries << endl;
	      
	      result = _fijk -> nEntries;
	    }
	}
    }
  
  return result;
}

int TtTagConsistencyEfficiency::readFijk( vector<string> mcFileName, vector<double> weight, string dataType, int minValue )
{
  int result = 0;

  string buf;
  struct Fijk * _fijk = NULL;
  if ( dataType == "signal" ) _fijk = &Fijk_sig;
  else if ( dataType == "background" ) _fijk = &Fijk_bg;
  else
    {
      cout << "Unknown data type, exiting..." << endl;
      result = -1;
    }
  if ( _fijk )
    {
      _fijk -> nEntries = 0;
      _fijk -> sum = 0;
      
      int i,j,k;
      int value = 0;
      const char * _format = "F_%d_%d_%d = %d";
      
      vector<string>::const_iterator _filename;
      for ( _filename = mcFileName . begin(); _filename != mcFileName . end(); _filename++  )
	{
	  ifstream _mcFile( _filename -> c_str() );
	    
	  if ( _mcFile . is_open() )
	    {
	      
	      while ( getline( _mcFile, buf ) > 0 )
		{
		  sscanf( buf . c_str(), _format, &i, &j, &k, &value );
		  //TEMPORARY!!!
		  //if ( value >= minValue )
		  if ( value >= minValue && (i+j+k) >= 4)
		    {
		      int _index = -1;
		      for ( int _i = 0; (_i < _fijk -> nEntries) && _index < 0; _i++ )
			{
			  if ( (_fijk->i[_i] == i) &&(_fijk->j[_i] == j) && (_fijk->k[_i] == k) ) _index = _i;
			}
		      if ( _index < 0 )
			{
			  _fijk -> i . push_back( i );
			  _fijk -> j . push_back( j );
			  _fijk -> k . push_back( k );
			  _fijk -> value[i][j][k] = (int)( (double)value * weight[_filename-mcFileName.begin()] );
			  _fijk -> nEntries++;
			}
		      else
			{
			  _fijk -> value[i][j][k] += (int)( (double)value * weight[_filename-mcFileName.begin()] );
			}
		      _fijk -> sum += value;
		    }
		}
	      
	      cout << "Number of Fijk read: " << _fijk -> nEntries << endl;
	      
	      result = _fijk -> nEntries;
	    }
	}
    }
  
  return result;
}

int TtTagConsistencyEfficiency::readFijk( vector<string> mcFileName, string dataType, int minValue )
{
  int result = 0;


  string buf;
  struct Fijk * _fijk = NULL;
  if ( dataType == "signal" ) _fijk = &Fijk_sig;
  else if ( dataType == "background" ) _fijk = &Fijk_bg;
  else
    {
      cout << "Unknown data type, exiting..." << endl;
      result = -1;
    }
  if ( _fijk )
    {
      _fijk -> nEntries = 0;
      _fijk -> sum = 0;
      
      int i,j,k;
      int value = 0;
      const char * _format = "F_%d_%d_%d = %d";
      
      vector<string>::const_iterator _filename;
      for ( _filename = mcFileName . begin(); _filename != mcFileName . end(); _filename++  )
	{
	  ifstream _mcFile( _filename -> c_str() );
	    
	  if ( _mcFile . is_open() )
	    {
	      
	      while ( getline( _mcFile, buf ) > 0 )
		{
		  sscanf( buf . c_str(), _format, &i, &j, &k, &value );
		  // TEMPORARY!!!
		  //if ( value >= minValue )
		  if ( value >= minValue && (i+j+k) >= 4)
		    {
		      int _index = -1;
		      for ( int _i = 0; (_i < _fijk -> nEntries) && _index < 0; _i++ )
			{
			  if ( (_fijk->i[_i] == i) &&(_fijk->j[_i] == j) && (_fijk->k[_i] == k) ) _index = _i;
			}
		      if ( _index < 0 )
			{
			  _fijk -> i . push_back( i );
			  _fijk -> j . push_back( j );
			  _fijk -> k . push_back( k );
			  _fijk -> value[i][j][k] = value;
			  _fijk -> nEntries++;
			}
		      else
			{
			  _fijk -> value[i][j][k] += value;
			}
		      _fijk -> sum += value;
		    }
		}
	      
	      cout << "Number of Fijk read: " << _fijk -> nEntries << endl;
	      
	      result = _fijk -> nEntries;
	    }
	}
    }
  
  return result;
}

int TtTagConsistencyEfficiency::readFijk( const char * mcFileName, string dataType, int minValue )
{
  int result = 0;

  ifstream _mcFile( mcFileName );
  string buf;
  struct Fijk * _fijk = NULL;
  if ( dataType == "signal" ) _fijk = &Fijk_sig;
  else if ( dataType == "background" ) _fijk = &Fijk_bg;
  else
    {
      cout << "Unknown data type, exiting..." << endl;
      result = -1;
    }
  if ( _fijk )
    {
      _fijk -> nEntries = 0;
      _fijk -> sum = 0;
      
      int i,j,k;
      int value = 0;
      const char * _format = "F_%d_%d_%d = %d";
      
      if ( _mcFile . is_open() )
	{
	  
	  while ( getline( _mcFile, buf ) > 0 )
	    {
	      sscanf( buf . c_str(), _format, &i, &j, &k, &value );
	      // TEMPORARY, REVISE!!!
	      //if ( value >= minValue )
	      if ( value >= minValue && (i+j+k) >= 4)
		{
		  _fijk -> i . push_back( i );
		  _fijk -> j . push_back( j );
		  _fijk -> k . push_back( k );
		  _fijk -> value[i][j][k] = value;
		  _fijk -> nEntries++;
		  _fijk -> sum += value;
		}
	    }
	  
	  cout << "Number of Fijk read: " << _fijk -> nEntries << endl;
	  
	  result = _fijk -> nEntries;
	}
      
    }

  return result;
}


int TtTagConsistencyEfficiency::sumFijk( const char * mcFileName )
{
  int result = -1;

  ifstream _mcFile( mcFileName );
  string buf;
  int sum = 0;
  int sum_underCut = 0;
  int sum_overCut = 0;
  int _fijk = 0;
  int i,j,k;
  const char * _format = "F_%d_%d_%d = %d";

  if ( _mcFile . is_open() )
    {
      
      while ( getline( _mcFile, buf ) > 0 )
	{
	  sscanf( buf . c_str(), _format, &i, &j, &k, &_fijk );
	  sum += _fijk;
	  if ( i + j + k < 4 ) sum_underCut += _fijk;
	  else sum_overCut += _fijk;
	}

      cout << "sum_underCut = " << sum_underCut << endl;
      cout << "sum_overCut = " << sum_overCut << endl;

      result = sum;
    }

  return result;
}

int TtTagConsistencyEfficiency::getN( int index )
{
  int result = -1;

  if ( inputFileRead )
    {
      result = N[index];
    }
  else
    {
      cout << "TtTagConsistencyEfficiency: input file has not been read. Exiting..." << endl;
    }

  return result;
}

int TtTagConsistencyEfficiency::getNb( void )
{
  int result = -1;

  if ( inputFileRead )
    {
      result = Nb;
    }
  else
    {
      cout << "TtTagConsistencyEfficiency: input file has not been read. Exiting..." << endl;
    }

  return result;
}

int TtTagConsistencyEfficiency::getNc( void )
{
  int result = -1;

  if ( inputFileRead )
    {
      result = Nc;
    }
  else
    {
      cout << "TtTagConsistencyEfficiency: input file has not been read. Exiting..." << endl;
    }

  return result;
}

int TtTagConsistencyEfficiency::getNl( void )
{
  int result = -1;

  if ( inputFileRead )
    {
      result = Nl;
    }
  else
    {
      cout << "TtTagConsistencyEfficiency: input file has not been read. Exiting..." << endl;
    }

  return result;
}

int TtTagConsistencyEfficiency::getNbtag( void )
{
  int result = -1;

  if ( inputFileRead )
    {
      result = Nbtag;
    }
  else
    {
      cout << "TtTagConsistencyEfficiency: input file has not been read. Exiting..." << endl;
    }

  return result;
}

int TtTagConsistencyEfficiency::getNctag( void )
{
  int result = -1;

  if ( inputFileRead )
    {
      result = Nctag;
    }
  else
    {
      cout << "TtTagConsistencyEfficiency: input file has not been read. Exiting..." << endl;
    }

  return result;
}

int TtTagConsistencyEfficiency::getNltag( void )
{
  int result = -1;

  if ( inputFileRead )
    {
      result = Nltag;
    }
  else
    {
      cout << "TtTagConsistencyEfficiency: input file has not been read. Exiting..." << endl;
    }

  return result;
}

int TtTagConsistencyEfficiency::readTable( const char * tabFileName, double discr )
{

  int result = 0;
  ifstream _tabFile( tabFileName );
  string buf;
  double aDiscr = -1.0;
  const char * _format = "%lf %d %d %d %d %d %d %d %d %d %d %d%*[^\n]";

  for (int i = 0; i < 5; i++ ) N[i] = 0;
  
  if ( _tabFile . is_open() )
    {

      do  // read line by line
	{
	  getline( _tabFile, buf );
	  sscanf( buf . c_str(), _format, &aDiscr, &N[0], &N[1], &N[2], &N[3], &N[4], &Nb, &Nbtag, &Nc, &Nctag, &Nl, &Nltag );
	} while( fabs(aDiscr - discr) > 0.0001 );

      inputFileRead = true;
      result = 0;
    }
  else
    {
      cout << "TtTagConsistencyEfficiency::readTable: cannot open input file" << endl;
      result = -1;
    }

  //cout << "discr = " << aDiscr << ", N_0 = " << N[0] << ", N_1 = " << N[1] << endl;

  return result;  
}

int TtTagConsistencyEfficiency::readTable( vector<string> tabFileName, double discr )
{

  int result = 0;
  //ifstream _tabFile( tabFileName );
  string buf;
  const char * _format = "%lf %d %d %d %d %d %d %d %d %d %d %d%*[^\n]";

  int _N[5];
  int _Nb, _Nc, _Nl, _Nbtag, _Nctag, _Nltag;

  for (int i = 0; i < 5; i++ )
    {
      N[i] = 0;
      _N[i] = 0;
    }
  Nb = 0;
  Nc = 0;
  Nl = 0;
  Nbtag = 0;
  Nctag = 0;
  Nltag = 0;
  _Nb = 0;
  _Nc = 0;
  _Nl = 0;
  _Nbtag = 0;
  _Nctag = 0;
  _Nltag = 0;

  
  vector<string>::const_iterator _filename;
  for ( _filename = tabFileName . begin(); _filename != tabFileName . end(); _filename++  )
    {
      ifstream _tabFile( _filename -> c_str() );
      double aDiscr = -1.0;
      
      if ( _tabFile . is_open() )
	{
	  
	  do  // read line by line
	    {
	      getline( _tabFile, buf );
	      sscanf( buf . c_str(), _format, &aDiscr, &_N[0], &_N[1], &_N[2], &_N[3], &_N[4], &_Nb, &_Nbtag, &_Nc, &_Nctag, &_Nl, &_Nltag );
	    } while( fabs(aDiscr - discr) > 0.0001 );
	  
	  for (int i = 0; i < 5; i++ ) N[i] += _N[i];
	  Nb += _Nb;
	  Nc += _Nc;
	  Nl += _Nl;
	  Nbtag += _Nbtag;
	  Nctag += _Nctag;
	  Nltag += _Nltag;      
	  
	  if ( fabs(aDiscr - discr) < 0.0001 )
	    {
	      inputFileRead = true;
	      result = 0;
	    }
	  else
	    {
	      cout << "TtTagConsistencyEfficiency::readTable: tagging discriminant doesn't match" << endl;
	      result = -1;
	    }
	}
      else
	{
	  cout << "TtTagConsistencyEfficiency::readTable: cannot open input file" << endl;
	  result = -1;
	}
      
    }
  
  //cout << "discr = " << aDiscr << ", N_0 = " << N[0] << ", N_1 = " << N[1] << endl;

  return result;  
}

int TtTagConsistencyEfficiency::updateTable( string tabListFileName, double discr, double weight, tableContent & _tab )
{

  int result = 0;
  string buf;
  const char * _format = "%lf %d %d %d %d %d %d %d %d %d %d %d%*[^\n]";

  int _N[5];
  int _Nb, _Nc, _Nl, _Nbtag, _Nctag, _Nltag;

  for (int i = 0; i < 5; i++ )
    {
      _N[i] = 0;
    }
  _Nb = 0;
  _Nc = 0;
  _Nl = 0;
  _Nbtag = 0;
  _Nctag = 0;
  _Nltag = 0;
  
  vector<string> tabFileName;
  char filename[1024];
  ifstream inFile( tabListFileName . c_str(), ios::in );
  if (!inFile)
    {
      cout << " Unable to open list file: " << tabListFileName << endl;
      return kFALSE;
    }
  else
    {
      cout << "List file opened successfully: " << tabListFileName << endl;
    }
  while (inFile >> filename)
    {
      string fullFileName = "";
      fullFileName . append( filename );
      tabFileName . push_back( fullFileName );
    }
  inFile.close();

  vector<string>::const_iterator _filename;
  for ( _filename = tabFileName . begin(); _filename != tabFileName . end(); _filename++  )
    {
      ifstream _tabFile( _filename -> c_str() );
      double aDiscr = -1.0;
      
      if ( _tabFile . is_open() )
	{
	  do  // read line by line
	    {
	      getline( _tabFile, buf );
	      sscanf( buf . c_str(), _format, &aDiscr, &_N[0], &_N[1], &_N[2], &_N[3], &_N[4], &_Nb, &_Nbtag, &_Nc, &_Nctag, &_Nl, &_Nltag );
	    } while( fabs(aDiscr - discr) > 0.0001 );
	  
	  if ( fabs(aDiscr - discr) < 0.0001 )
	    {
	      for (int i = 0; i < 5; i++ ) _tab . N[i] += (int)( (double)_N[i] * weight );
	      _tab . Nb += (int)( (double)_Nb * weight );
	      _tab . Nc += (int)( (double)_Nc * weight );
	      _tab . Nl += (int)( (double)_Nl * weight );
	      _tab . Nbtag += (int)( (double)_Nbtag * weight );
	      _tab . Nctag += (int)( (double)_Nctag * weight );
	      _tab . Nltag += (int)( (double)_Nltag * weight );      
	  
	      result = 0;
	    }
	  else
	    {
	      cout << "TtTagConsistencyEfficiency::readTable: tagging discriminant doesn't match" << endl;
	      result = -1;
	    }
	}
      else
	{
	  cout << "TtTagConsistencyEfficiency::readTable: cannot open input file" << endl;
	  result = -1;
	}
      
    }
  
  cout << "discr = " << discr << ", N_1 = " << _tab.N[1] << ", N_2 = " << _tab.N[2] << ", N_3= " << _tab.N[3] << endl;

  return result;  
}

int TtTagConsistencyEfficiency::readTable( vector<string> tabFileName, double discr, vector<double> weight, tableContent & _tab )
{

  int result = 0;
  //ifstream _tabFile( tabFileName );
  string buf;
  const char * _format = "%lf %d %d %d %d %d %d %d %d %d %d %d%*[^\n]";

  int _N[5];
  int _Nb, _Nc, _Nl, _Nbtag, _Nctag, _Nltag;

  for (int i = 0; i < 5; i++ )
    {
      _tab . N[i] = 0;
      _N[i] = 0;
    }
  _tab . Nb = 0;
  _tab . Nc = 0;
  _tab . Nl = 0;
  _tab . Nbtag = 0;
  _tab . Nctag = 0;
  _tab . Nltag = 0;
  _Nb = 0;
  _Nc = 0;
  _Nl = 0;
  _Nbtag = 0;
  _Nctag = 0;
  _Nltag = 0;

  
  vector<string>::const_iterator _filename;
  for ( _filename = tabFileName . begin(); _filename != tabFileName . end(); _filename++  )
    {
      ifstream _tabFile( _filename -> c_str() );
      double aDiscr = -1.0;
      
      if ( _tabFile . is_open() )
	{
	  do  // read line by line
	    {
	      getline( _tabFile, buf );
	      sscanf( buf . c_str(), _format, &aDiscr, &_N[0], &_N[1], &_N[2], &_N[3], &_N[4], &_Nb, &_Nbtag, &_Nc, &_Nctag, &_Nl, &_Nltag );
	    } while( fabs(aDiscr - discr) > 0.0001 );
	  
	  if ( fabs(aDiscr - discr) < 0.0001 )
	    {
	      for (int i = 0; i < 5; i++ ) _tab . N[i] += (int)( (double)_N[i] * weight[_filename-tabFileName.begin()] );
	      _tab . Nb += (int)( (double)_Nb * weight[_filename-tabFileName.begin()] );
	      _tab . Nc += (int)( (double)_Nc * weight[_filename-tabFileName.begin()] );
	      _tab . Nl += (int)( (double)_Nl * weight[_filename-tabFileName.begin()] );
	      _tab . Nbtag += (int)( (double)_Nbtag * weight[_filename-tabFileName.begin()] );
	      _tab . Nctag += (int)( (double)_Nctag * weight[_filename-tabFileName.begin()] );
	      _tab . Nltag += (int)( (double)_Nltag * weight[_filename-tabFileName.begin()] );      
	  
	      result = 0;
	    }
	  else
	    {
	      cout << "TtTagConsistencyEfficiency::readTable: tagging discriminant doesn't match" << endl;
	      result = -1;
	    }
	}
      else
	{
	  cout << "TtTagConsistencyEfficiency::readTable: cannot open input file" << endl;
	  result = -1;
	}
      
    }
  
  //cout << "discr = " << aDiscr << ", N_0 = " << N[0] << ", N_1 = " << N[1] << endl;

  return result;  
}

int TtTagConsistencyEfficiency::updateTable( vector<string> tabFileName, double discr, double weight, tableLine & _tab )
{

  int result = 0;
  //ifstream _tabFile( tabFileName );
  string buf;
  const char * _format = "%lf %d %d %d %d %d %d %d %d %d %d %d %d %d%*[^\n]";

  int _N[5];
  int _Nb, _Nc, _Nl, _Nx, _Nbtag, _Nctag, _Nltag, _Nxtag;

  for (int i = 0; i < 5; i++ )
    {
      _N[i] = 0;
    }
  _Nb = 0;
  _Nc = 0;
  _Nl = 0;
  _Nx = 0;
  _Nbtag = 0;
  _Nctag = 0;
  _Nltag = 0;
  _Nxtag = 0;

  
  vector<string>::const_iterator _filename;
  for ( _filename = tabFileName . begin(); _filename != tabFileName . end(); _filename++  )
    {
      ifstream _tabFile( _filename -> c_str() );
      double aDiscr = -1.0;
      
      if ( _tabFile . is_open() )
	{
	  do  // read line by line
	    {
	      getline( _tabFile, buf );
	      sscanf( buf . c_str(), _format, &aDiscr, &_N[0], &_N[1], &_N[2], &_N[3], &_N[4], &_Nb, &_Nbtag, &_Nc, &_Nctag, &_Nl, &_Nltag, &_Nx, &_Nxtag );
	    } while( fabs(aDiscr - discr) > 0.0001 );
	  
	  if ( fabs(aDiscr - discr) < 0.0001 )
	    {
	      for (int i = 0; i < 5; i++ ) _tab . N[i] += (double)_N[i] * weight;
	      _tab . discr = aDiscr;
	      _tab . Nb += (double)_Nb * weight;
	      _tab . Nc += (double)_Nc * weight;
	      _tab . Nl += (double)_Nl * weight;
	      _tab . Nx += (double)_Nx * weight;
	      _tab . Nbtag += (double)_Nbtag * weight;
	      _tab . Nctag += (double)_Nctag * weight;
	      _tab . Nltag += (double)_Nltag * weight;      
	      _tab . Nxtag += (double)_Nxtag * weight;      
	  
	      result = 0;
	    }
	  else
	    {
	      cout << "TtTagConsistencyEfficiency::readTable: tagging discriminant doesn't match" << endl;
	      result = -1;
	    }
	}
      else
	{
	  cout << "TtTagConsistencyEfficiency::readTable: cannot open input file" << endl;
	  result = -1;
	}
      
    }
  
  //cout << "discr = " << aDiscr << ", N_0 = " << N[0] << ", N_1 = " << N[1] << endl;

  return result;  
}

int TtTagConsistencyEfficiency::updateTable( string tabFileName, double discr, double weight, tableLine & _tab )
{

  int result = 0;
  //ifstream _tabFile( tabFileName );
  string buf;
  const char * _format = "%lf %d %d %d %d %d %d %d %d %d %d %d %d %d%*[^\n]";

  int _N[5];
  int _Nb, _Nc, _Nl, _Nx, _Nbtag, _Nctag, _Nltag, _Nxtag;

  for (int i = 0; i < 5; i++ )
    {
      _N[i] = 0;
    }
  _Nb = 0;
  _Nc = 0;
  _Nl = 0;
  _Nx = 0;
  _Nbtag = 0;
  _Nctag = 0;
  _Nltag = 0;
  _Nxtag = 0;


  string * _filename = &tabFileName;
  
      ifstream _tabFile( _filename -> c_str() );
      double aDiscr = -1.0;
      
      if ( _tabFile . is_open() )
	{
	  do  // read line by line
	    {
	      getline( _tabFile, buf );
	      sscanf( buf . c_str(), _format, &aDiscr, &_N[0], &_N[1], &_N[2], &_N[3], &_N[4], &_Nb, &_Nbtag, &_Nc, &_Nctag, &_Nl, &_Nltag, &_Nx, &_Nxtag );
	    } while( fabs(aDiscr - discr) > 0.0001 );
	  
	  if ( fabs(aDiscr - discr) < 0.0001 )
	    {
	      for (int i = 0; i < 5; i++ ) _tab . N[i] += (double)_N[i] * weight;
	      _tab . discr = aDiscr;
	      _tab . Nb += (double)_Nb * weight;
	      _tab . Nc += (double)_Nc * weight;
	      _tab . Nl += (double)_Nl * weight;
	      _tab . Nx += (double)_Nx * weight;
	      _tab . Nbtag += (double)_Nbtag * weight;
	      _tab . Nctag += (double)_Nctag * weight;
	      _tab . Nltag += (double)_Nltag * weight;      
	      _tab . Nxtag += (double)_Nxtag * weight;      
	  
	      result = 0;
	    }
	  else
	    {
	      cout << "TtTagConsistencyEfficiency::readTable: tagging discriminant doesn't match" << endl;
	      result = -1;
	    }
	}
      else
	{
	  cout << "TtTagConsistencyEfficiency::readTable: cannot open input file" << endl;
	  result = -1;
	}
      
  //cout << "discr = " << aDiscr << ", N_0 = " << N[0] << ", N_1 = " << N[1] << endl;

  return result;  
}




int TtTagConsistencyEfficiency::updateTableFromList( string tabListFileName, double discr, double weight, tableLine & _tab )
{

  int result = 0;
  //ifstream _tabFile( tabFileName );
  string buf;
  const char * _format = "%lf %d %d %d %d %d %d %d %d %d %d %d %d %d%*[^\n]";

  int _N[5];
  int _Nb, _Nc, _Nl, _Nx, _Nbtag, _Nctag, _Nltag, _Nxtag;

  for (int i = 0; i < 5; i++ )
    {
      _N[i] = 0;
    }
  _Nb = 0;
  _Nc = 0;
  _Nl = 0;
  _Nx = 0;
  _Nbtag = 0;
  _Nctag = 0;
  _Nltag = 0;
  _Nxtag = 0;


  //  string * _filename = &tabFileName;
  
  vector<string> tabFileName;
  char filename[1024];
  ifstream inFile( tabListFileName . c_str(), ios::in );
  if (!inFile)
    {
      cout << " Unable to open list file: " << tabListFileName << endl;
      return kFALSE;
    }
  else
    {
      cout << "List file opened successfully: " << tabListFileName << endl;
    }
  while (inFile >> filename)
    {
      string fullFileName = "";
      fullFileName . append( filename );
      tabFileName . push_back( fullFileName );
    }
  inFile.close();
  
  vector<string>::const_iterator _filename;
  for ( _filename = tabFileName . begin(); _filename != tabFileName . end(); _filename++  )
    {
      ifstream _tabFile( _filename -> c_str() );
      double aDiscr = -1.0;
      
      if ( _tabFile . is_open() )
	{
	  do  // read line by line
	    {
	      getline( _tabFile, buf );
	      sscanf( buf . c_str(), _format, &aDiscr, &_N[0], &_N[1], &_N[2], &_N[3], &_N[4], &_Nb, &_Nbtag, &_Nc, &_Nctag, &_Nl, &_Nltag, &_Nx, &_Nxtag );
	    } while( fabs(aDiscr - discr) > 0.0001 );
	  
	  if ( fabs(aDiscr - discr) < 0.0001 )
	    {
	      for (int i = 0; i < 5; i++ ) _tab . N[i] += (double)_N[i] * weight;
	      _tab . discr = aDiscr;
	      _tab . Nb += (double)_Nb * weight;
	      _tab . Nc += (double)_Nc * weight;
	      _tab . Nl += (double)_Nl * weight;
	      _tab . Nx += (double)_Nx * weight;
	      _tab . Nbtag += (double)_Nbtag * weight;
	      _tab . Nctag += (double)_Nctag * weight;
	      _tab . Nltag += (double)_Nltag * weight;      
	      _tab . Nxtag += (double)_Nxtag * weight;      
	  
	      result = 0;
	    }
	  else
	    {
	      cout << "TtTagConsistencyEfficiency::readTable: tagging discriminant doesn't match" << endl;
	      result = -1;
	    }
	}
      else
	{
	  cout << "TtTagConsistencyEfficiency::readTable: cannot open input file" << endl;
	  result = -1;
	}
    }
  cout << "discr = " << discr << ", N_1 = " << _tab.N[1] << ", N_2 = " << _tab.N[2] << ", N_3= " << _tab.N[3] << endl;

  return result;  
}

RooFitResult TtTagConsistencyEfficiency::minTest2( int n0, int n1, int n2, int n3, int n4, const char * option,  double eb, double ec, double el, double xsec_ttbar, double bgf, double eb_min, double eb_max, double ec_min, double ec_max, double el_min, double el_max, double xsec_ttbar_min, double xsec_ttbar_max, double bgf_min, double bgf_max )
{

  RooRealVar nEvt1( "nEvt1", "observed number of tags", n1, 0.9*n1, 1.1*n1 ); // number of events with 1 tag
  RooRealVar nEvt2( "nEvt2", "observed number of tags", n2, 0.9*n2, 1.1*n2 );
  RooRealVar nEvt3( "nEvt3", "observed number of tags", n3, 0.9*n3, 1.1*n3 );

  RooRealVar lumi( "lumi", "Integrated luminocity", 100, 0, 1000 );     // 1/pb
  RooRealVar xsec( "xsec", "ttbar cross section", xsec_ttbar, xsec_ttbar_min, xsec_ttbar_max );      // pb
  //RooRealVar  eff( "eff",  "Pre-tag efficiency",    0.10625, 0.0, 1.0 );
  RooRealVar  eff( "eff",  "Pre-tag efficiency",    0.169, 0.0, 1.0 );

  RooRealVar bgFrac( "bgFrac",  "Fraction of background", bgf, bgf_min, bgf_max );

  RooRealVar epsb( "epsb", "b tag efficiency", eb, max(eb_min,0.0), min(eb_max,1.0) );  // b tagging efficiency
  RooRealVar epsc( "epsc", "c tag efficiency", ec, max(ec_min,0.0), min(ec_max,1.0) );  // c tagging efficiency
  RooRealVar epsl( "epsl", "l tag efficiency", el, max(el_min,0.0), min(el_max,1.0) );  // light (and gluon) tagging efficiency - taken from outside !!!

  RooArgList lamargs[5];

  for ( int m = 0; m < 5; m++ )
    {
      lamargs[m] . add(lumi);
      lamargs[m] . add(xsec);
      lamargs[m] . add(eff);
    }

  //
  // ===> Signal coefficients
  //
  map<int, map<int, map<int, RooRealVar *> > > _varFijk;
  map<int, map<int, map<int, map<int, RooFormulaVar *> > > > _Fijkn;
  TtTagConsistencyEfficiency::BFormula _f;
  ostringstream lam_sumFsig[5];
  for ( int l = 1; l < 4; l++ ) lam_sumFsig[l] << "(";
  for ( unsigned int ii = 0; ii < Fijk_sig.i.size(); ii++ )
    {
      int _i = Fijk_sig.i[ii];
      int _j = Fijk_sig.j[ii];
      int _k = Fijk_sig.k[ii];
      _f = getFormula( _i, _j, _k );
      _varFijk[_i][_j][_k] = new RooRealVar( _f . name . c_str(), _f . title . c_str(), (double)Fijk_sig.value[_i][_j][_k]/(double)Fijk_sig.sum, 0.0, 1.0 );
      for ( int l = 1; l < 4; l++ )
	{
	  ostringstream _name;
	  _name . str("");
	  _name << _f . name . c_str() << "_" << l; 
	  if ( lam_sumFsig[l] . str() . size() > 1 ) lam_sumFsig[l] << "+";
	  lam_sumFsig[l] << _name . str();
	  if ( _i == 0 && _j != 0 && _k != 0 ) _Fijkn[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk[_i][_j][_k], epsc, epsl ) );
	  if ( _i != 0 && _j == 0 && _k != 0 ) _Fijkn[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk[_i][_j][_k], epsb, epsl ) );
	  if ( _i != 0 && _j != 0 && _k == 0 ) _Fijkn[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk[_i][_j][_k], epsb, epsc ) );
	  if ( _i == 0 && _j == 0 && _k != 0 ) _Fijkn[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk[_i][_j][_k], epsl ) );
	  if ( _i == 0 && _j != 0 && _k == 0 ) _Fijkn[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk[_i][_j][_k], epsc ) );
	  if ( _i != 0 && _j == 0 && _k == 0 ) _Fijkn[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk[_i][_j][_k], epsb ) );
	  else _Fijkn[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk[_i][_j][_k], epsb, epsc, epsl ) );
	  lamargs[l] . add( *_Fijkn[_i][_j][_k][l] );
	  //_Fijkn[_i][_j][_k][l] -> Print();
	}
      _varFijk[_i][_j][_k] -> setConstant( kTRUE );
      //_varFijk[_i][_j][_k] -> Print();
    }
  for ( int l = 1; l < 4; l++ ) lam_sumFsig[l] << ")";

  for ( int m = 0; m < 5; m++ ) lamargs[m] . add(bgFrac);
  
  //
  // =====> Background coefficients
  //
  map<int, map<int, map<int, RooRealVar *> > > _varFijk_bg;
  map<int, map<int, map<int, map<int, RooFormulaVar *> > > > _Fijkn_bg;
  ostringstream lam_sumFbg[5];
  for ( int l = 1; l < 4; l++ ) lam_sumFbg[l] << "(";
  for ( unsigned int ii = 0; ii < Fijk_bg.i.size(); ii++ )
    {
      int _i = Fijk_bg.i[ii];
      int _j = Fijk_bg.j[ii];
      int _k = Fijk_bg.k[ii];
      _f = getFormula( _i, _j, _k, "_bg" );
      _varFijk_bg[_i][_j][_k] = new RooRealVar( _f . name . c_str(), _f . title . c_str(), (double)Fijk_bg.value[_i][_j][_k]/(double)Fijk_bg.sum, 0.0, 1.0 );
      for ( int l = 1; l < 4; l++ )
	{
	  ostringstream _name;
	  _name . str("");
	  _name << _f . name . c_str() << "_" << l; 
	  if ( lam_sumFbg[l] . str() . size() > 1 ) lam_sumFbg[l] << "+";
	  lam_sumFbg[l] << _name . str();
	  if ( _i == 0 && _j != 0 && _k != 0 ) _Fijkn_bg[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk_bg[_i][_j][_k], epsc, epsl ) );
	  if ( _i != 0 && _j == 0 && _k != 0 ) _Fijkn_bg[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk_bg[_i][_j][_k], epsb, epsl ) );
	  if ( _i != 0 && _j != 0 && _k == 0 ) _Fijkn_bg[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk_bg[_i][_j][_k], epsb, epsc ) );
	  if ( _i == 0 && _j == 0 && _k != 0 ) _Fijkn_bg[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk_bg[_i][_j][_k], epsl ) );
	  if ( _i == 0 && _j != 0 && _k == 0 ) _Fijkn_bg[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk_bg[_i][_j][_k], epsc ) );
	  if ( _i != 0 && _j == 0 && _k == 0 ) _Fijkn_bg[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk_bg[_i][_j][_k], epsb ) );
	  else _Fijkn_bg[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk_bg[_i][_j][_k], epsb, epsc, epsl ) );
	  lamargs[l] . add( *_Fijkn_bg[_i][_j][_k][l] );
	  //_Fijkn_bg[_i][_j][_k][l] -> Print();
	}
      _varFijk_bg[_i][_j][_k] -> setConstant( kTRUE );
      //_varFijk_bg[_i][_j][_k] -> Print();
    }
  for ( int l = 1; l < 4; l++ ) lam_sumFbg[l] << ")";

  //
  // ===========================>

  string lam1_str = "lumi*xsec*eff*(" + lam_sumFsig[1] . str() + "+bgFrac*" + lam_sumFbg[1] . str() + ")";
  string lam2_str = "lumi*xsec*eff*(" + lam_sumFsig[2] . str() + "+bgFrac*" + lam_sumFbg[2] . str() + ")";
  string lam3_str = "lumi*xsec*eff*(" + lam_sumFsig[3] . str() + "+bgFrac*" + lam_sumFbg[3] . str() + ")";

  RooFormulaVar lam1( "lam1", lam1_str . c_str(), lamargs[1] );
  RooFormulaVar lam2( "lam2", lam2_str . c_str(), lamargs[2] );
  RooFormulaVar lam3( "lam3", lam3_str . c_str(), lamargs[3] );

  //lamargs[1] . Print();
  //lam1 . Print();

  lumi . setConstant(kTRUE);
  eff  . setConstant(kTRUE);
  bgFrac . setConstant(kTRUE);
  epsl . setConstant(kTRUE);

  nEvt1 . setConstant(kTRUE);
  nEvt2 . setConstant(kTRUE);
  nEvt3 . setConstant(kTRUE);

  //
  // ===> minimization
  //
  RooFormulaVar poisson1( "poisson1", "exp(-(nEvt1-lam1)*(nEvt1-lam1)/nEvt1/nEvt1)", RooArgList(nEvt1,lam1) );
  RooFormulaVar poisson2( "poisson2", "exp(-(nEvt2-lam2)*(nEvt2-lam2)/nEvt2/nEvt2)", RooArgList(nEvt2,lam2) );
  RooFormulaVar poisson3( "poisson3", "exp(-(nEvt3-lam3)*(nEvt3-lam3)/nEvt3/nEvt3)", RooArgList(nEvt3,lam3) );

  RooFormulaVar LL( "LL", "-log(poisson1*poisson2*poisson3)", RooArgList(poisson1,poisson2,poisson3) );  

  RooMinuit m(LL);

  m . fit( "msh" );

  RooFormulaVar p1( "p1", "TMath::Gaus( (nEvt1-lam1)/sqrt(lam1), 0.0, 1.0 )", RooArgList(nEvt1,lam1) );
  RooFormulaVar p2( "p2", "TMath::Gaus( (nEvt2-lam2)/sqrt(lam2), 0.0, 1.0 )", RooArgList(nEvt2,lam2) );
  RooFormulaVar p3( "p3", "TMath::Gaus( (nEvt3-lam3)/sqrt(lam3), 0.0, 1.0 )", RooArgList(nEvt3,lam3) );

  RooFormulaVar LL2( "LL2", "-log(p1*p2*p3)", RooArgList(p1,p2,p3) );  

  RooMinuit m2(LL2);

  RooFitResult * result = m2 . fit( option );

  //contourBC = m2 . contour( epsb, epsc, 1, 2, 3 );
  //contourBX = m2 . contour( epsb, xsec, 1, 2, 3 );

  return *result;
}




RooFitResult TtTagConsistencyEfficiency::fit( TtTagConsistencyFitConfig & config )
{
  //cout << "DEBUG: n1,2,3 = " << config.n1 << "   " << config . n2 << "   " << config.n3 << endl; 


  /*
  int _temp = Fijk_sig . value[1][0][3];
  Fijk_sig . value[1][0][3] += _temp * 0.25;
  Fijk_sig . sum += _temp * 0.25;

  _temp = Fijk_sig . value[1][0][4];
  Fijk_sig . value[1][0][4] += _temp * 0.25;
  Fijk_sig . sum += _temp * 0.25;

  _temp = Fijk_sig . value[2][0][3];
  Fijk_sig . value[2][0][3] += _temp * 0.2;
  Fijk_sig . sum += _temp * 0.22;
  */




  RooRealVar nEvt1( "nEvt1", "observed number of tags", config . n1, 0.9*config . n1, 1.1*config . n1 ); // number of events with 1 tag
  RooRealVar nEvt2( "nEvt2", "observed number of tags", config . n2, 0.9*config . n2, 1.1*config . n2 );
  RooRealVar nEvt3( "nEvt3", "observed number of tags", config . n3, 0.9*config . n3, 1.1*config . n3 );

  RooRealVar lumi( "lumi", "Integrated luminocity", config . lumi, 0, 1000 );     // 1/pb
  RooRealVar xsec( "xsec", "ttbar cross section", config . xsec_ttbar, config . xsec_ttbar_min, config . xsec_ttbar_max );      // pb
  RooRealVar  eff( "eff",  "Pre-tag efficiency", config . efficiency, 0.0, 1.0 );

  RooRealVar bgFrac0( "bgFrac0",  "Fraction of background in n0 bin", config . bgf0, config . bgf0_min, config . bgf0_max );
  RooRealVar bgFrac1( "bgFrac1",  "Fraction of background in n1 bin", config . bgf1, config . bgf1_min, config . bgf1_max );
  RooRealVar bgFrac2( "bgFrac2",  "Fraction of background in n2 bin", config . bgf2, config . bgf2_min, config . bgf2_max );
  RooRealVar bgFrac3( "bgFrac3",  "Fraction of background in n3 bin", config . bgf3, config . bgf3_min, config . bgf3_max );
  RooRealVar bgFrac4( "bgFrac4",  "Fraction of background in n4 bin", config . bgf4, config . bgf4_min, config . bgf4_max );

  RooRealVar epsb( "epsb", "b tag efficiency", config . eb,  max(config . eb_min,0.0), min(config . eb_max,1.0) );  // b tagging efficiency
  RooRealVar epsc( "epsc", "c tag efficiency", config . ec,  max(config . ec_min,0.0), min(config . ec_max,1.0) );  // c tagging efficiency
  RooRealVar epsl( "epsl", "l tag efficiency", config . el,  max(config . el_min,0.0), min(config . el_max,1.0) );  // light (and gluon) tagging efficiency - from outside !

  RooArgList lamargs[5];

  for ( int m = 0; m < 5; m++ )
    {
      lamargs[m] . add(lumi);
      lamargs[m] . add(xsec);
      lamargs[m] . add(eff);
    }

  //
  // ===> Signal coefficients
  //
  map<int, map<int, map<int, RooRealVar *> > > _varFijk;
  map<int, map<int, map<int, map<int, RooFormulaVar *> > > > _Fijkn;
  TtTagConsistencyEfficiency::BFormula _f;
  ostringstream lam_sumFsig[5];
  for ( int l = 1; l < 4; l++ ) lam_sumFsig[l] << "(";
  for ( unsigned int ii = 0; ii < Fijk_sig.i.size(); ii++ )
    {
      int _i = Fijk_sig.i[ii];
      int _j = Fijk_sig.j[ii];
      int _k = Fijk_sig.k[ii];
      _f = getFormula( _i, _j, _k );
      _varFijk[_i][_j][_k] = new RooRealVar( _f . name . c_str(), _f . title . c_str(), (double)Fijk_sig.value[_i][_j][_k]/(double)Fijk_sig.sum, 0.0, 1.0 );
      for ( int l = 1; l < 4; l++ )
	{
	  ostringstream _name;
	  _name . str("");
	  _name << _f . name . c_str() << "_" << l; 
	  if ( lam_sumFsig[l] . str() . size() > 1 ) lam_sumFsig[l] << "+";
	  lam_sumFsig[l] << _name . str();
	  if ( _i == 0 && _j != 0 && _k != 0 ) _Fijkn[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk[_i][_j][_k], epsc, epsl ) );
	  if ( _i != 0 && _j == 0 && _k != 0 ) _Fijkn[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk[_i][_j][_k], epsb, epsl ) );
	  if ( _i != 0 && _j != 0 && _k == 0 ) _Fijkn[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk[_i][_j][_k], epsb, epsc ) );
	  if ( _i == 0 && _j == 0 && _k != 0 ) _Fijkn[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk[_i][_j][_k], epsl ) );
	  if ( _i == 0 && _j != 0 && _k == 0 ) _Fijkn[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk[_i][_j][_k], epsc ) );
	  if ( _i != 0 && _j == 0 && _k == 0 ) _Fijkn[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk[_i][_j][_k], epsb ) );
	  else _Fijkn[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk[_i][_j][_k], epsb, epsc, epsl ) );
	  lamargs[l] . add( *_Fijkn[_i][_j][_k][l] );
	  //_Fijkn[_i][_j][_k][l] -> Print();
	}
      _varFijk[_i][_j][_k] -> setConstant( kTRUE );
      //_varFijk[_i][_j][_k] -> Print();
    }
  for ( int l = 1; l < 4; l++ ) lam_sumFsig[l] << ")";

  //for ( int m = 0; m < 5; m++ ) lamargs[m] . add(bgFrac);
  lamargs[0] . add(bgFrac0);
  lamargs[1] . add(bgFrac1);
  lamargs[2] . add(bgFrac2);
  lamargs[3] . add(bgFrac3);
  lamargs[4] . add(bgFrac4);
  
  //
  // =====> Background coefficients
  //
  map<int, map<int, map<int, RooRealVar *> > > _varFijk_bg;
  map<int, map<int, map<int, map<int, RooFormulaVar *> > > > _Fijkn_bg;
  ostringstream lam_sumFbg[5];
  for ( int l = 1; l < 4; l++ ) lam_sumFbg[l] << "(";
  for ( unsigned int ii = 0; ii < Fijk_bg.i.size(); ii++ )
    {
      int _i = Fijk_bg.i[ii];
      int _j = Fijk_bg.j[ii];
      int _k = Fijk_bg.k[ii];
      _f = getFormula( _i, _j, _k, "_bg" );
      _varFijk_bg[_i][_j][_k] = new RooRealVar( _f . name . c_str(), _f . title . c_str(), (double)Fijk_bg.value[_i][_j][_k]/(double)Fijk_bg.sum, 0.0, 1.0 );
      for ( int l = 1; l < 4; l++ )
	{
	  ostringstream _name;
	  _name . str("");
	  _name << _f . name . c_str() << "_" << l; 
	  if ( lam_sumFbg[l] . str() . size() > 1 ) lam_sumFbg[l] << "+";
	  lam_sumFbg[l] << _name . str();
	  if ( _i == 0 && _j != 0 && _k != 0 ) _Fijkn_bg[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk_bg[_i][_j][_k], epsc, epsl ) );
	  if ( _i != 0 && _j == 0 && _k != 0 ) _Fijkn_bg[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk_bg[_i][_j][_k], epsb, epsl ) );
	  if ( _i != 0 && _j != 0 && _k == 0 ) _Fijkn_bg[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk_bg[_i][_j][_k], epsb, epsc ) );
	  if ( _i == 0 && _j == 0 && _k != 0 ) _Fijkn_bg[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk_bg[_i][_j][_k], epsl ) );
	  if ( _i == 0 && _j != 0 && _k == 0 ) _Fijkn_bg[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk_bg[_i][_j][_k], epsc ) );
	  if ( _i != 0 && _j == 0 && _k == 0 ) _Fijkn_bg[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk_bg[_i][_j][_k], epsb ) );
	  else _Fijkn_bg[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk_bg[_i][_j][_k], epsb, epsc, epsl ) );
	  lamargs[l] . add( *_Fijkn_bg[_i][_j][_k][l] );
	  //_Fijkn_bg[_i][_j][_k][l] -> Print();
	}
      _varFijk_bg[_i][_j][_k] -> setConstant( kTRUE );
      //_varFijk_bg[_i][_j][_k] -> Print();
    }
  for ( int l = 1; l < 4; l++ ) lam_sumFbg[l] << ")";

  //
  // ===========================>

  string lam1_str = "lumi*xsec*eff*(" + lam_sumFsig[1] . str() + "+bgFrac1*" + lam_sumFbg[1] . str() + ")";
  string lam2_str = "lumi*xsec*eff*(" + lam_sumFsig[2] . str() + "+bgFrac2*" + lam_sumFbg[2] . str() + ")";
  string lam3_str = "lumi*xsec*eff*(" + lam_sumFsig[3] . str() + "+bgFrac3*" + lam_sumFbg[3] . str() + ")";

  RooFormulaVar lam1( "lam1", lam1_str . c_str(), lamargs[1] );
  RooFormulaVar lam2( "lam2", lam2_str . c_str(), lamargs[2] );
  RooFormulaVar lam3( "lam3", lam3_str . c_str(), lamargs[3] );

  //lamargs[1] . Print();
  //lam1 . Print();

  lumi . setConstant(kTRUE);
  eff  . setConstant(kTRUE);
  bgFrac1 . setConstant(kTRUE);
  bgFrac2 . setConstant(kTRUE);
  bgFrac3 . setConstant(kTRUE);
  //epsc . setConstant(kTRUE);
  epsl . setConstant(kTRUE);

  //xsec . setConstant(kTRUE);

  nEvt1 . setConstant(kTRUE);
  nEvt2 . setConstant(kTRUE);
  nEvt3 . setConstant(kTRUE);

  //
  // ===> minimization
  //
  RooFormulaVar poisson1( "poisson1", "exp(-(nEvt1-lam1)*(nEvt1-lam1)/nEvt1/nEvt1)", RooArgList(nEvt1,lam1) );
  RooFormulaVar poisson2( "poisson2", "exp(-(nEvt2-lam2)*(nEvt2-lam2)/nEvt2/nEvt2)", RooArgList(nEvt2,lam2) );
  RooFormulaVar poisson3( "poisson3", "exp(-(nEvt3-lam3)*(nEvt3-lam3)/nEvt3/nEvt3)", RooArgList(nEvt3,lam3) );

  RooFormulaVar LL( "LL", "-log(poisson1*poisson2*poisson3)", RooArgList(poisson1,poisson2,poisson3) );  

  RooMinuit m(LL);

  m . fit( "msh" );

  RooFormulaVar p1( "p1", "TMath::Gaus( (nEvt1-lam1)/sqrt(lam1), 0.0, 1.0 )", RooArgList(nEvt1,lam1) );
  RooFormulaVar p2( "p2", "TMath::Gaus( (nEvt2-lam2)/sqrt(lam2), 0.0, 1.0 )", RooArgList(nEvt2,lam2) );
  RooFormulaVar p3( "p3", "TMath::Gaus( (nEvt3-lam3)/sqrt(lam3), 0.0, 1.0 )", RooArgList(nEvt3,lam3) );

  RooFormulaVar LL2( "LL2", "-log(p1*p2*p3)", RooArgList(p1,p2,p3) );  

  RooMinuit m2(LL2);

  RooFitResult * result = m2 . fit( config . option );

  //contourBC = m2 . contour( epsb, epsc, 1, 2, 3 );
  //contourBX = m2 . contour( epsb, xsec, 1, 2, 3 );

  return *result;
}


RooFitResult TtTagConsistencyEfficiency::fit_signal( TtTagConsistencyFitConfig & config, TH2F * _contour )
{
  //cout << "DEBUG: n1,2,3 = " << config.n1 << "   " << config . n2 << "   " << config.n3 << endl; 


  RooRealVar nEvt1( "nEvt1", "observed number of tags", config . n1, 0.9*config . n1, 1.1*config . n1 ); // number of events with 1 tag
  RooRealVar nEvt2( "nEvt2", "observed number of tags", config . n2, 0.9*config . n2, 1.1*config . n2 );
  RooRealVar nEvt3( "nEvt3", "observed number of tags", config . n3, 0.9*config . n3, 1.1*config . n3 );

  RooRealVar lumi( "lumi", "Integrated luminocity", config . lumi, 0, 1000 );     // 1/pb
  RooRealVar xsec( "xsec", "ttbar cross section", config . xsec_ttbar, config . xsec_ttbar_min, config . xsec_ttbar_max );      // pb
  RooRealVar  eff( "eff",  "Pre-tag efficiency", config . efficiency, 0.0, 1.0 );

  RooRealVar epsb( "epsb", "b tag efficiency", config . eb,  max(config . eb_min,0.0), min(config . eb_max,1.0) );  // b tagging efficiency
  RooRealVar epsc( "epsc", "c tag efficiency", config . ec,  max(config . ec_min,0.0), min(config . ec_max,1.0) );  // c tagging efficiency
  RooRealVar epsl( "epsl", "l tag efficiency", config . el,  max(config . el_min,0.0), min(config . el_max,1.0) );  // light (and gluon) tagging efficiency - from outside !

  RooArgList lamargs[5];

  for ( int m = 0; m < 5; m++ )
    {
      lamargs[m] . add(lumi);
      lamargs[m] . add(xsec);
      lamargs[m] . add(eff);
    }

  //
  // ===> Signal coefficients
  //
  map<int, map<int, map<int, RooRealVar *> > > _varFijk;
  map<int, map<int, map<int, map<int, RooFormulaVar *> > > > _Fijkn;
  TtTagConsistencyEfficiency::BFormula _f;
  ostringstream lam_sumFsig[5];
  for ( int l = 1; l < 4; l++ ) lam_sumFsig[l] << "(";
  for ( unsigned int ii = 0; ii < Fijk_sig.i.size(); ii++ )
    {
      int _i = Fijk_sig.i[ii];
      int _j = Fijk_sig.j[ii];
      int _k = Fijk_sig.k[ii];
      _f = getFormula( _i, _j, _k );
      _varFijk[_i][_j][_k] = new RooRealVar( _f . name . c_str(), _f . title . c_str(), (double)Fijk_sig.value[_i][_j][_k]/(double)Fijk_sig.sum, 0.0, 1.0 );
      for ( int l = 1; l < 4; l++ )
	{
	  ostringstream _name;
	  _name . str("");
	  _name << _f . name . c_str() << "_" << l; 
	  if ( lam_sumFsig[l] . str() . size() > 1 ) lam_sumFsig[l] << "+";
	  lam_sumFsig[l] << _name . str();
	  if ( _i == 0 && _j != 0 && _k != 0 ) _Fijkn[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk[_i][_j][_k], epsc, epsl ) );
	  if ( _i != 0 && _j == 0 && _k != 0 ) _Fijkn[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk[_i][_j][_k], epsb, epsl ) );
	  if ( _i != 0 && _j != 0 && _k == 0 ) _Fijkn[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk[_i][_j][_k], epsb, epsc ) );
	  if ( _i == 0 && _j == 0 && _k != 0 ) _Fijkn[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk[_i][_j][_k], epsl ) );
	  if ( _i == 0 && _j != 0 && _k == 0 ) _Fijkn[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk[_i][_j][_k], epsc ) );
	  if ( _i != 0 && _j == 0 && _k == 0 ) _Fijkn[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk[_i][_j][_k], epsb ) );
	  else _Fijkn[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk[_i][_j][_k], epsb, epsc, epsl ) );
	  lamargs[l] . add( *_Fijkn[_i][_j][_k][l] );
	  //_Fijkn[_i][_j][_k][l] -> Print();
	}
      _varFijk[_i][_j][_k] -> setConstant( kTRUE );
      //_varFijk[_i][_j][_k] -> Print();
    }
  for ( int l = 1; l < 4; l++ ) lam_sumFsig[l] << ")";
  
  string lam1_str = "lumi*xsec*eff*(" + lam_sumFsig[1] . str()  + ")";
  string lam2_str = "lumi*xsec*eff*(" + lam_sumFsig[2] . str()  + ")";
  string lam3_str = "lumi*xsec*eff*(" + lam_sumFsig[3] . str()  + ")";

  RooFormulaVar lam1( "lam1", lam1_str . c_str(), lamargs[1] );
  RooFormulaVar lam2( "lam2", lam2_str . c_str(), lamargs[2] );
  RooFormulaVar lam3( "lam3", lam3_str . c_str(), lamargs[3] );

  //lamargs[1] . Print();
  //lam1 . Print();

  lumi . setConstant(kTRUE);
  eff  . setConstant(kTRUE);
  epsl . setConstant(kTRUE);

  //xsec . setConstant(kTRUE);

  nEvt1 . setConstant(kTRUE);
  nEvt2 . setConstant(kTRUE);
  nEvt3 . setConstant(kTRUE);

  //
  // ===> minimization
  //
  RooFormulaVar poisson1( "poisson1", "exp(-(nEvt1-lam1)*(nEvt1-lam1)/nEvt1/nEvt1)", RooArgList(nEvt1,lam1) );
  RooFormulaVar poisson2( "poisson2", "exp(-(nEvt2-lam2)*(nEvt2-lam2)/nEvt2/nEvt2)", RooArgList(nEvt2,lam2) );
  RooFormulaVar poisson3( "poisson3", "exp(-(nEvt3-lam3)*(nEvt3-lam3)/nEvt3/nEvt3)", RooArgList(nEvt3,lam3) );

  RooFormulaVar LL( "LL", "-log(poisson1*poisson2*poisson3)", RooArgList(poisson1,poisson2,poisson3) );  

  RooMinuit m(LL);

  m . fit( "msh" );

  RooFormulaVar p1( "p1", "TMath::Gaus( (nEvt1-lam1)/sqrt(lam1), 0.0, 1.0 )", RooArgList(nEvt1,lam1) );
  RooFormulaVar p2( "p2", "TMath::Gaus( (nEvt2-lam2)/sqrt(lam2), 0.0, 1.0 )", RooArgList(nEvt2,lam2) );
  RooFormulaVar p3( "p3", "TMath::Gaus( (nEvt3-lam3)/sqrt(lam3), 0.0, 1.0 )", RooArgList(nEvt3,lam3) );

  RooFormulaVar LL2( "LL2", "-log(p1*p2*p3)", RooArgList(p1,p2,p3) );  

  RooMinuit m2(LL2);

  RooFitResult * result = m2 . fit( config . option );

  //contourBC = m2 . contour( epsb, epsc, 1, 2, 3 );
  //contourBX = m2 . contour( epsb, xsec, 1, 2, 3 );

  if ( _contour != NULL )
    {
      *_contour = *(m2 . contour( epsb, epsc, 1, 2, 3 ));
    }

  return *result;
}





RooFitResult TtTagConsistencyEfficiency::fit_optimized( TtTagConsistencyFitConfig & config )
{

  RooRealVar nEvt1( "nEvt1", "observed number of tags", config . n1, 0.9*config . n1, 1.1*config . n1 ); // number of events with 1 tag
  RooRealVar nEvt2( "nEvt2", "observed number of tags", config . n2, 0.9*config . n2, 1.1*config . n2 );
  RooRealVar nEvt3( "nEvt3", "observed number of tags", config . n3, 0.9*config . n3, 1.1*config . n3 );

  RooRealVar lumi( "lumi", "Integrated luminocity", config . lumi, 0, 1000 );     // 1/pb
  RooRealVar xsec( "xsec", "ttbar cross section", config . xsec_ttbar, config . xsec_ttbar_min, config . xsec_ttbar_max );      // pb
  RooRealVar  eff( "eff",  "Pre-tag efficiency", config . efficiency, 0.0, 1.0 );

  RooRealVar bgFrac0( "bgFrac0",  "Fraction of background in n0 bin", config . bgf0, config . bgf0_min, config . bgf0_max );
  RooRealVar bgFrac1( "bgFrac1",  "Fraction of background in n1 bin", config . bgf1, config . bgf1_min, config . bgf1_max );
  RooRealVar bgFrac2( "bgFrac2",  "Fraction of background in n2 bin", config . bgf2, config . bgf2_min, config . bgf2_max );
  RooRealVar bgFrac3( "bgFrac3",  "Fraction of background in n3 bin", config . bgf3, config . bgf3_min, config . bgf3_max );
  RooRealVar bgFrac4( "bgFrac4",  "Fraction of background in n4 bin", config . bgf4, config . bgf4_min, config . bgf4_max );

  RooRealVar epsb( "epsb", "b tag efficiency", config . eb,  max(config . eb_min,0.0), min(config . eb_max,1.0) );  // b tagging efficiency
  RooRealVar epsc( "epsc", "c tag efficiency", config . ec,  max(config . ec_min,0.0), min(config . ec_max,1.0) );  // c tagging efficiency
  RooRealVar epsl( "epsl", "l tag efficiency", config . el,  max(config . el_min,0.0), min(config . el_max,1.0) );  // light (and gluon) tagging efficiency - from outside !

  RooArgList lamargs[5];

  for ( int m = 0; m < 5; m++ )
    {
      lamargs[m] . add(lumi);
      lamargs[m] . add(xsec);
      lamargs[m] . add(eff);
    }

  //
  // ===> Signal coefficients
  //
  map<int, map<int, map<int, RooRealVar *> > > _varFijk;
  map<int, map<int, map<int, map<int, RooFormulaVar *> > > > _Fijkn;
  TtTagConsistencyEfficiency::BFormula _f;
  ostringstream lam_sumFsig[5];
  for ( int l = 1; l < 4; l++ ) lam_sumFsig[l] << "(";
  for ( unsigned int ii = 0; ii < Fijk_sig.i.size(); ii++ )
    {
      int _i = Fijk_sig.i[ii];
      int _j = Fijk_sig.j[ii];
      int _k = Fijk_sig.k[ii];
      _f = getFormula( _i, _j, _k );
      _varFijk[_i][_j][_k] = new RooRealVar( _f . name . c_str(), _f . title . c_str(), (double)Fijk_sig.value[_i][_j][_k]/(double)Fijk_sig.sum, 0.0, 1.0 );
      for ( int l = 1; l < 4; l++ )
	{
	  ostringstream _name;
	  _name . str("");
	  _name << _f . name . c_str() << "_" << l; 
	  if ( lam_sumFsig[l] . str() . size() > 1 ) lam_sumFsig[l] << "+";
	  lam_sumFsig[l] << _name . str();
	  if ( _i == 0 && _j != 0 && _k != 0 ) _Fijkn[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk[_i][_j][_k], epsc, epsl ) );
	  if ( _i != 0 && _j == 0 && _k != 0 ) _Fijkn[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk[_i][_j][_k], epsb, epsl ) );
	  if ( _i != 0 && _j != 0 && _k == 0 ) _Fijkn[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk[_i][_j][_k], epsb, epsc ) );
	  if ( _i == 0 && _j == 0 && _k != 0 ) _Fijkn[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk[_i][_j][_k], epsl ) );
	  if ( _i == 0 && _j != 0 && _k == 0 ) _Fijkn[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk[_i][_j][_k], epsc ) );
	  if ( _i != 0 && _j == 0 && _k == 0 ) _Fijkn[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk[_i][_j][_k], epsb ) );
	  else _Fijkn[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk[_i][_j][_k], epsb, epsc, epsl ) );
	  lamargs[l] . add( *_Fijkn[_i][_j][_k][l] );
	  //_Fijkn[_i][_j][_k][l] -> Print();
	}
      _varFijk[_i][_j][_k] -> setConstant( kTRUE );
      //_varFijk[_i][_j][_k] -> Print();
    }
  for ( int l = 1; l < 4; l++ ) lam_sumFsig[l] << ")";

  //for ( int m = 0; m < 5; m++ ) lamargs[m] . add(bgFrac);
  lamargs[0] . add(bgFrac0);
  lamargs[1] . add(bgFrac1);
  lamargs[2] . add(bgFrac2);
  lamargs[3] . add(bgFrac3);
  lamargs[4] . add(bgFrac4);
  
  //
  // =====> Background coefficients
  //
  map<int, map<int, map<int, RooRealVar *> > > _varFijk_bg;
  map<int, map<int, map<int, map<int, RooFormulaVar *> > > > _Fijkn_bg;
  ostringstream lam_sumFbg[5];
  for ( int l = 1; l < 4; l++ ) lam_sumFbg[l] << "(";
  for ( unsigned int ii = 0; ii < Fijk_bg.i.size(); ii++ )
    {
      int _i = Fijk_bg.i[ii];
      int _j = Fijk_bg.j[ii];
      int _k = Fijk_bg.k[ii];
      _f = getFormula( _i, _j, _k, "_bg" );
      _varFijk_bg[_i][_j][_k] = new RooRealVar( _f . name . c_str(), _f . title . c_str(), (double)Fijk_bg.value[_i][_j][_k]/(double)Fijk_bg.sum, 0.0, 1.0 );
      for ( int l = 1; l < 4; l++ )
	{
	  ostringstream _name;
	  _name . str("");
	  _name << _f . name . c_str() << "_" << l; 
	  if ( lam_sumFbg[l] . str() . size() > 1 ) lam_sumFbg[l] << "+";
	  lam_sumFbg[l] << _name . str();
	  if ( _i == 0 && _j != 0 && _k != 0 ) _Fijkn_bg[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk_bg[_i][_j][_k], epsc, epsl ) );
	  if ( _i != 0 && _j == 0 && _k != 0 ) _Fijkn_bg[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk_bg[_i][_j][_k], epsb, epsl ) );
	  if ( _i != 0 && _j != 0 && _k == 0 ) _Fijkn_bg[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk_bg[_i][_j][_k], epsb, epsc ) );
	  if ( _i == 0 && _j == 0 && _k != 0 ) _Fijkn_bg[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk_bg[_i][_j][_k], epsl ) );
	  if ( _i == 0 && _j != 0 && _k == 0 ) _Fijkn_bg[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk_bg[_i][_j][_k], epsc ) );
	  if ( _i != 0 && _j == 0 && _k == 0 ) _Fijkn_bg[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk_bg[_i][_j][_k], epsb ) );
	  else _Fijkn_bg[_i][_j][_k][l] = new RooFormulaVar( _name . str() . c_str(), _f . formula[l] . c_str(), RooArgList( *_varFijk_bg[_i][_j][_k], epsb, epsc, epsl ) );
	  lamargs[l] . add( *_Fijkn_bg[_i][_j][_k][l] );
	  //_Fijkn_bg[_i][_j][_k][l] -> Print();
	}
      _varFijk_bg[_i][_j][_k] -> setConstant( kTRUE );
      //_varFijk_bg[_i][_j][_k] -> Print();
    }
  for ( int l = 1; l < 4; l++ ) lam_sumFbg[l] << ")";

  //
  // ===========================>

  string lam1_str = "lumi*xsec*eff*(" + lam_sumFsig[1] . str() + "+bgFrac1*" + lam_sumFbg[1] . str() + ")";
  string lam2_str = "lumi*xsec*eff*(" + lam_sumFsig[2] . str() + "+bgFrac2*" + lam_sumFbg[2] . str() + ")";
  string lam3_str = "lumi*xsec*eff*(" + lam_sumFsig[3] . str() + "+bgFrac3*" + lam_sumFbg[3] . str() + ")";

  RooFormulaVar lam1( "lam1", lam1_str . c_str(), lamargs[1] );
  RooFormulaVar lam2( "lam2", lam2_str . c_str(), lamargs[2] );
  RooFormulaVar lam3( "lam3", lam3_str . c_str(), lamargs[3] );

  //lamargs[1] . Print();
  //lam1 . Print();

  lumi . setConstant(kTRUE);
  eff  . setConstant(kTRUE);
  bgFrac1 . setConstant(kTRUE);
  bgFrac2 . setConstant(kTRUE);
  bgFrac3 . setConstant(kTRUE);
  epsc . setConstant(kTRUE);
  epsl . setConstant(kTRUE);

  //xsec . setConstant(kTRUE);

  nEvt1 . setConstant(kTRUE);
  nEvt2 . setConstant(kTRUE);
  nEvt3 . setConstant(kTRUE);

  //
  // ===> minimization
  //
  RooFormulaVar poisson1( "poisson1", "exp(-(nEvt1-lam1)*(nEvt1-lam1)/nEvt1/nEvt1)", RooArgList(nEvt1,lam1) );
  RooFormulaVar poisson2( "poisson2", "exp(-(nEvt2-lam2)*(nEvt2-lam2)/nEvt2/nEvt2)", RooArgList(nEvt2,lam2) );
  RooFormulaVar poisson3( "poisson3", "exp(-(nEvt3-lam3)*(nEvt3-lam3)/nEvt3/nEvt3)", RooArgList(nEvt3,lam3) );

  RooFormulaVar LL( "LL", "-log(poisson1*poisson2*poisson3)", RooArgList(poisson1,poisson2,poisson3) );  

  RooMinuit m(LL);

  m . fit( "msh" );

  RooFormulaVar p1( "p1", "TMath::Gaus( (nEvt1-lam1)/sqrt(lam1), 0.0, 1.0 )", RooArgList(nEvt1,lam1) );
  RooFormulaVar p2( "p2", "TMath::Gaus( (nEvt2-lam2)/sqrt(lam2), 0.0, 1.0 )", RooArgList(nEvt2,lam2) );
  RooFormulaVar p3( "p3", "TMath::Gaus( (nEvt3-lam3)/sqrt(lam3), 0.0, 1.0 )", RooArgList(nEvt3,lam3) );

  RooFormulaVar LL2( "LL2", "-log(p1*p2*p3)", RooArgList(p1,p2,p3) );  

  RooMinuit m2(LL2);

  RooFitResult * result = m2 . fit( config . option );

  //contourBC = m2 . contour( epsb, epsc, 1, 2, 3 );
  //contourBX = m2 . contour( epsb, xsec, 1, 2, 3 );

  return *result;
}


