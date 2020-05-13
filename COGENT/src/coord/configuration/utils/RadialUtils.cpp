#include <tuple>
#include <string>
#include <list>
#include <fstream>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include "RadialUtils.H"
#include "Directions.H"
#include "Simulation.H"

#include "NamespaceHeader.H"

bool fillDataBRadial(HermiteSpline& a_spline, std::string& a_str)
{
  /** The method creates a Hermite cubic spline based on the data read from file,
   *  and returns it implicitly as "a_spline" reference.
   *  File name is a_str, and the file is processed by readDataBRadial().
   *  All correct data readings are performed in the readDataBRadial().
   */
  
  int N;  /// Size of list of pairs read from file
  /// Yes, we create a useless list for all processes but procID()==0, but it is empty anyway
  std::list< std::pair<double, double> > lst;
  if (!procID()) {
    readDataBRadial(lst, a_str);
    N = lst.size();
  }
  
  /// Broadcast N to all to all the processes, so that they can allocate memory now
  MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
  double* ax = new double[N];
  double* ay = new double[N];
  
  /// Read the file and fill list of pairs for 0-th process
  if (!procID()) {
    /// Fill the arrays ax and ay
    std::list< std::pair<double, double> >::iterator it;
    int counter = 0;
    for(it=lst.begin(); it!=lst.end(); ++it){
      ax[counter] = it->first;
      ay[counter] = it->second;
      counter++;
    }
    /// Check if the array ax is sorted
    for(int i=0;i<N-1;i++){
      if (ax[i+1]<ax[i]) {
	delete[] ax;
	delete[] ay;
	MayDay::Error("ERROR: fillDataBRadial() failed! Data is not sorted.");}
    }
  }
  
  /// Send these arrays to all the processes
  MPI_Bcast(ax, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(ay, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  /// Define the spline and delete arrays
  a_spline.buildSpline(ax, ay, N, false);
  
  delete[] ax;
  delete[] ay;
  return true;
}

bool readDataBRadial(std::list< std::pair<double, double> >& a_lst,  std::string& a_str)
{
  /** The method opens the file with a_str filename and reads the data out from it.
   *  It is assumed that the data has a format: "double separator double",
   *  where separator is any number of ' ' or '\t' characters.
   *  The method is also protected from some bad inpuit,
   *  as any number of ' ' and '\t' characters are allowed around any number that is read.
   *  The method calls MAyDay::Error() if something does not work nicely.
   */
  
  std::fstream file_read(a_str.c_str());
  if (!file_read.is_open() ) {
    string tmp_str = "ERROR: readDataBRadial() failed! File \"";
    tmp_str += a_str;
    tmp_str += "\" cannot be opened.";
    MayDay::Error(tmp_str.c_str());
  }
  
  std::string line;

  /// Loop over all the file and collect all the lines
  while(getline(file_read, line)) {
    /// Sweep through the line to get both values
    /// Two words times two points (start, finish). Initialized all to -1 in order to easily check if the input is read correctly
    int position[4] = {-1, -1, -1, -1};
    int i0 = 0;   /// Start searching from zero position
    /// Loop over two words we are looking for
    for(int n_words=0; n_words<2; n_words++){
      bool word_started = false;
      int i;
      for(i=i0; i<line.size(); i++){
        if (!word_started && line[i]!=' ' && line[i]!='\t'){
          position[2*n_words] = i;
          word_started = true;
          continue;
        }
        if ( word_started && (line[i]==' ' || line[i]=='\t') ){
          position[2*n_words+1] = i;
          i0 = i+1;
          break;
        }
      }
      if (i==line.size()) {break;}
    }
    if (position[2]==-1) {MayDay::Error("ERROR: readDataBRadial() failed! Pairs of data are corrupt.");}
    if (position[3]==-1) {position[3] = line.size();}   /// Because the line ended exactly at the last char of the second number

    /// Create a pair and add it to the list
    /// Validity of the input is checked by std::stod
    std::string tmp_str;
    tmp_str = line.substr(position[0], position[1]-position[0]);
    std::pair<double, double> tmp_pair;
    
    try {tmp_pair.first = std::stod(tmp_str);}
    catch (const std::invalid_argument&) {
      MayDay::Error("ERROR: readDataBRadial() failed! File contains corrupt data.");
      throw;
    }
    
    tmp_str = line.substr(position[2], position[3]-position[2]);
    try {tmp_pair.second = std::stod(tmp_str);}
    catch (const std::invalid_argument&) {
      MayDay::Error("ERROR: readDataBRadial() failed! File contains corrupt data.");
      throw;
    }
    
    a_lst.push_back(tmp_pair);
  }
  
  if (a_lst.size()==0) {MayDay::Error("ERROR: readDataBRadial() failed! File contains no data.");}
  return true;
}


/// Empty constructor
HermiteSpline::HermiteSpline()
{
  m_initialized = false;
  m_a = nullptr;
  m_b = nullptr;
  m_c = nullptr;
  m_h = nullptr;
  m_arr_x = nullptr;
  m_arr_y = nullptr;
  m_N = -1;
}

/// Constructor with data
HermiteSpline::HermiteSpline(double* a_ax, double* a_ay, int a_N_, bool a_derivatives, double a_a_0_, double a_a_N_, double a_err)
{
  /** m_arr_x[i] is the SORTED array of x-nodes
   *  m_arr_y[i] is the array of y-values
   *  "derivatives" is a flag if derivatives at the boundaries are given
   *  a_N_ is the total number of points, INCLUDING end points. There are a_N_-1 domains.
   */
  m_initialized = false;
  m_a = nullptr;
  m_b = nullptr;
  m_c = nullptr;
  m_h = nullptr;
  m_arr_x = nullptr;
  m_arr_y = nullptr;
  m_N = -1;
  this->buildSpline(a_ax, a_ay, a_N_, a_derivatives, a_a_0_, a_a_N_, a_err);

}

bool HermiteSpline::buildSpline(double* a_ax, double* a_ay, int a_N_, bool a_derivatives, double a_a_0_, double a_a_N_, double a_err)
{
  /** a_ax[i] is the SORTED array of x-nodes
   *  a_ay[i] is the array of y-values
   *  "derivatives" is a flag if derivatives at the boundaries are given
   *  a_N_ is the total number of points, INCLUDING end points. There are a_N_-1 domains.
   */
  
  /// Check for sorted property
  for(int i=0; i<a_N_-1 ; i++){
    if(a_ax[i+1] < a_ax[i]) {std::cerr << "ERROR: HermiteSpline::buildSpline() failed! Array of the coordinates is not sorted!\n"; return false;}
  }
  
  if(m_N!= a_N_-1){       /// Check if we already have spline created
    if(!m_initialized){  /// If not, create new one
      m_N = a_N_-1;
      m_a = new double[m_N];
      m_b = new double[m_N];
      m_c = new double[m_N];
      m_h = new double[m_N];
      m_arr_x = new double[m_N+1];
      m_arr_y = new double[m_N+1];
    }
    else{           /// If yes, delete the previous one and create a new one
      delete[] m_a;
      delete[] m_b;
      delete[] m_c;
      delete[] m_h;
      delete[] m_arr_x;
      delete[] m_arr_y;
      
      m_N = a_N_-1;
      m_a = new double[m_N];
      m_b = new double[m_N];
      m_c = new double[m_N];
      m_h = new double[m_N];
      m_arr_x = new double[m_N+1];
      m_arr_y = new double[m_N+1];
    }
  }   /// Otherwise, the spline of the proper length is already created, so we just update it
  
  /// Dummy initialization in case the object is not constructed correctly
  /// Also proper filling of arr_x[] and arr_y[]
  for(int i=0;i<m_N;i++){
    m_a[i] = 2.3e205;
    m_b[i] = 2.3e205;
    m_c[i] = 2.3e205;
    m_h[i] = 2.3e205;
    m_arr_x[i] = a_ax[i];
    m_arr_y[i] = a_ay[i];
  }
  m_arr_x[m_N] = a_ax[m_N];
  m_arr_y[m_N] = a_ay[m_N];
  
  if(!a_derivatives){
    if(a_N_<4) {
      std::cerr << "ERROR: HermiteSpline::buildSpline() failed! Too few points provided, or derivatives are missing.\n";
      return false;
    }
    /// Derivative at x[0]
    double p_abc[3];
    cubicFit(p_abc, m_arr_x[0], m_arr_x, m_arr_y);
    m_a_0 = p_abc[0];
    /// Derivative at x[N]
    double m_arr_x4[4], m_arr_y4[4];
    for(int i=0;i<4;i++) {m_arr_x4[i] = m_arr_x[a_N_-1-i]; m_arr_y4[i] = m_arr_y[a_N_-1-i];}
    cubicFit(p_abc, m_arr_x[m_N], m_arr_x4, m_arr_y4);
    m_a_N = p_abc[0];
  }
  else {m_a_0 = a_a_0_; m_a_N = a_a_N_;}
  
  /// Main routine
  m_ERR = a_err;
  m_a[0] = m_a_0;
  for(int i=0;i<m_N;i++) {
    m_h[i] = m_arr_x[i+1] - m_arr_x[i];
    if(m_h[i]<m_ERR) {
      std::cerr << "ERROR: HermiteSpline::buildSpline() failed! Either x[i] array is not sorted, or the points are too close to each other.\n";
      return false;
    }
  }
  
  /// N=1 special case
  if(m_N==1){
    m_c[0] = ( (m_a_0+m_a_N)*m_h[0] - 2*(m_arr_y[1]-m_arr_y[0]) ) / m_h[0]/m_h[0]/m_h[0];
    m_b[0] = ( 3*(m_arr_y[1]-m_arr_y[0]) - m_h[0]*(2*m_a_0+m_a_N) ) / m_h[0]/m_h[0];
    return true;
  }
  /// N=2 special case
  if(m_N==2){
    m_a[1] = ( 3*m_h[0]*(m_arr_y[2]-m_arr_y[1])/m_h[1] + 3*m_h[1]*(m_arr_y[1]-m_arr_y[0])/m_h[0] - m_a[0]*m_h[1] - m_a_N*m_h[0] ) / (m_h[0]+m_h[1])/2;
    m_b[0] = 3*(m_arr_y[1]-m_arr_y[0])/m_h[0]/m_h[0] - (2*m_a[0]+m_a[1])/m_h[0];
    m_b[1] = 3*(m_arr_y[2]-m_arr_y[1])/m_h[1]/m_h[1] - (2*m_a[1]+m_a_N)/m_h[1];
    m_c[0] = (m_b[1]-m_b[0])/3/m_h[0];
    m_c[1] = (m_a_N+m_a[1])/m_h[1]/m_h[1] + 2*(m_arr_y[1]-m_arr_y[2])/m_h[1]/m_h[1]/m_h[1];
    return true;
  }
  
  /// N>2 case
  double* f = new double[m_N-1];    /// array for RHS
  double* s = new double[m_N-1];    /// array for the center diagonal
  double* t = new double[m_N-2];    /// array for the upper diagonal
  double* v = new double[m_N-2];    /// array for the lower diagonal
  
  /// Constructing the matrix
  for(int i=0;i<m_N-2;i++){
    s[i] = 2*(m_h[i]+m_h[i+1]);
    t[i] = m_h[i];
    v[i] = m_h[i+2];
    f[i] = 3*(m_h[i]*(m_arr_y[i+2]-m_arr_y[i+1])/m_h[i+1] + m_h[i+1]*(m_arr_y[i+1]-m_arr_y[i])/m_h[i]);
  }
  f[0] = f[0] - m_a[0]*m_h[1];
  s[m_N-2] = 2*(m_h[m_N-2]+m_h[m_N-1]);
  f[m_N-2] = 3*(m_h[m_N-2]*(m_arr_y[m_N]-m_arr_y[m_N-1])/m_h[m_N-1] + m_h[m_N-1]*(m_arr_y[m_N-1]-m_arr_y[m_N-2])/m_h[m_N-2]) - m_a_N*m_h[m_N-2];
  /// Matrix inversion
  threeDiagSolver(m_a, s, t, v, f, m_N-1);
  /// Position it in a properly
  for(int i=m_N-1;i>0;i--) {m_a[i] = m_a[i-1];}
  m_a[0] = m_a_0;
  /// Set all a, b, c
  for(int i=0;i<m_N-1;i++){
    m_b[i] = 3*(m_arr_y[i+1]-m_arr_y[i])/m_h[i]/m_h[i] - (2*m_a[i] + m_a[i+1])/m_h[i];
    m_c[i] = (m_a[i+1]+m_a[i])/m_h[i]/m_h[i] + 2*(m_arr_y[i]-m_arr_y[i+1])/m_h[i]/m_h[i]/m_h[i];
  }
  m_c[m_N-1] = (m_a_N+m_a[m_N-1])/m_h[m_N-1]/m_h[m_N-1] + 2*(m_arr_y[m_N-1]-m_arr_y[m_N])/m_h[m_N-1]/m_h[m_N-1]/m_h[m_N-1];
  m_b[m_N-1] = 3*(m_arr_y[m_N]-m_arr_y[m_N-1])/m_h[m_N-1]/m_h[m_N-1] - (2*m_a[m_N-1] + m_a_N)/m_h[m_N-1];
  
  delete[] f;
  delete[] s;
  delete[] t;
  delete[] v;
  
  return true;
}

/// Destructor
HermiteSpline::~HermiteSpline()
{
  if(m_initialized){
    delete[] m_a;
    delete[] m_b;
    delete[] m_c;
    delete[] m_h;
    delete[] m_arr_x;
    delete[] m_arr_y;
  }
}

double HermiteSpline::getValAt(int a_ind, double a_x) const
{
  /// Returns delta spline value y-y[i] around point "ind" at the coordinate h[ind]>x>0
  if (a_ind<0 || a_ind>m_N-1) {std::cerr<<"ERROR: HermiteSpline::getValAt() failed! Index is out of range.\n"; return 1.3e209;}
  if (a_x<0 || a_x>m_h[a_ind]) {std::cerr<<"ERROR: HermiteSpline::getValAt() failed! Wrong range of \"x\".\n"; return 1.3e209;}
  return m_a[a_ind]*a_x + m_b[a_ind]*a_x*a_x + m_c[a_ind]*a_x*a_x*a_x;
}

void HermiteSpline::getValAt(double* a_res, double a_x) const
{
  /** The method implicitly returns the value of the spline interpolation at the position a_x.
   *  It also returns first and second derivatives at this point.
   *  a_res[0] = value, a_res[1] = first derivative, a_res[2] = second derivative.
   *  If the coordinate is out of the domain, the method linearly extrapolates data there
   */
  
  /// Check if a_x inside the domain
  if(a_x<m_arr_x[0] || a_x>m_arr_x[m_N]) {
    std::cerr << "Warning: HermiteSpline::getValAt()! Coordinate is out of range. Linear extrapolation is used.\n";
    std::cerr << "coord_min: " << m_arr_x[0] << "  coord_given: " << a_x << "  coord_max " << m_arr_x[m_N] << "\n";
    if (a_x>m_arr_x[m_N]) {
      a_res[2] = 0.0;
      a_res[1] = m_a_N;
      a_res[0] = m_arr_y[m_N] + (a_x - m_arr_x[m_N]) * a_res[1];
      return;
    }
    if (a_x<m_arr_x[0]) {
      a_res[2] = 0.0;
      a_res[1] = m_a_0;
      a_res[0] = m_arr_y[0] + (a_x - m_arr_x[0]) * a_res[1];
      return;
    }
  }
  
  /// Perform logarithmic search with initial guess of a homogeneous array arr_x
  double D = (a_x - m_arr_x[0]) / (m_arr_x[m_N] - m_arr_x[0]);
  int ind = int(m_N*D);
  int ind_max = m_N-1;
  int ind_min = 0;
  int counter = 0;
  while(true) {
    counter++;
    if (counter>100) {
      /// Iterative procedure failed, return some rubbish
      std::cout << "ERROR: HermiteSpline<data_type>::getValAt() failed!"
                << "Position or index are not found.\n"
                << ind << "\t\t" << ind_min << "\t\t" << ind_max << "\n"
                << m_arr_x[ind] << "\t\t" << a_x << "\t\t" << m_arr_x[ind+1] << "\n"
                << m_arr_x[0] << "\t\t" << m_arr_x[m_N] << "\n";
      a_res[0] = 1.3e209;
      a_res[1] = 1.3e209;
      a_res[2] = 1.3e209;
      return;
    }
    
    if (a_x<m_arr_x[ind]) {
      ind_max = ind;
      ind = (ind + ind_min)/2;
      continue;
    }
    if (a_x>m_arr_x[ind+1]) {
      ind_min = ind;
      ind = (ind + ind_max)/2;
      continue;
    }
    break;  /// Both conditions are satisfied, we leave the loop
  }
  D = a_x - m_arr_x[ind]; /// delta x from the nearest neighbor
  
  /// Return results
  a_res[0] = m_arr_y[ind] + m_a[ind]*D + m_b[ind]*D*D + m_c[ind]*D*D*D;
  a_res[1] = m_a[ind] + 2*m_b[ind]*D + 3*m_c[ind]*D*D;
  a_res[2] = 2*m_b[ind] + 6*m_c[ind]*D;
  return;
}

double cubicFit(double* a_p_abc, double a_x_f, double* a_x, double* a_y, double a_ERR)
{
  /** The function creates a third-order polynomial
   *  y = y[0] + a*(x-x[0]) + b*(x-x[0])**2 + c*(x-x[0])**3,
   *  based on 4 provided points for x and for y:
   *  y(x[0]) = a_y[0], y(x[1]) = a_y[1], y(x[2]) = a_y[2], y(x[3]) = a_y[3].
   *  It returns y_f = a_y[x_f] and also transforms the polynomial to x_f-centered one:
   *  y = y_f + a_p_abc[0]*(x-a_x_f) + a_p_abc[1]*(x-a_x_f)**2 + a_p_abc[2]*(x-a_x_f)**3,
   *  where a_p_abc[2] = c; a_p_abc[1] = b + 3*c*(a_x_f-a_x[0]); a_p_abc[0] = a + 2*b*(a_x_f-a_x[0]) + 3*c*(a_x_f-a_x[0])**2.
   *  The array p_abc[] is also implicitly returned.
   */
  double tmp;
  double D3 = (a_x[3]-a_x[2]) * (a_x[3]-a_x[1]) * (a_x[3]-a_x[0]);
  double D2 = (a_x[2]-a_x[1]) * (a_x[2]-a_x[0]) * (a_x[1]-a_x[0]);
  if(fabs(D3)<a_ERR || fabs(D2)<a_ERR) {
    std::cerr << "ERROR: cubicFit() failed! x-points are too close to each other!\n";
    /// Returns extremely large value in order to break the code
    return 1.0e200;
  }
  
  /// Doing "a" term by term
  tmp = (a_x[2]-a_x[0]) * (a_x[2]-a_x[0]) * (a_x[1]-a_x[0]) * (a_x[1]-a_x[0]) * (a_x[2]-a_x[1]);
  a_p_abc[0] = tmp*a_y[3];
  tmp = (a_x[3]-a_x[0]) * (a_x[3]-a_x[0]) * (a_x[1]-a_x[0]) * (a_x[1]-a_x[0]) * (a_x[3]-a_x[1]);
  a_p_abc[0] -= tmp*a_y[2];
  tmp = (a_x[3]-a_x[0]) * (a_x[3]-a_x[0]) * (a_x[2]-a_x[0]) * (a_x[2]-a_x[0]) * (a_x[3]-a_x[2]);
  a_p_abc[0] += tmp*a_y[1];
  tmp = (a_x[3]-a_x[2]) * (a_x[3]-a_x[1]) * (a_x[2]-a_x[1]);
  tmp *= ( (a_x[3]-a_x[0])*(a_x[2]-a_x[0]) + (a_x[3]-a_x[0])*(a_x[1]-a_x[0]) + (a_x[2]-a_x[0])*(a_x[1]-a_x[0]) );
  a_p_abc[0] -= tmp*a_y[0];
  
  /// Doing "b" term by term
  tmp = (a_x[2]-a_x[1]) * (a_x[2]-a_x[0]) * (a_x[1]-a_x[0]) * (a_x[2]+a_x[1]-2*a_x[0]);
  a_p_abc[1] = -tmp*a_y[3];
  tmp = (a_x[3]-a_x[1]) * (a_x[3]-a_x[0]) * (a_x[1]-a_x[0]) * (a_x[3]+a_x[1]-2*a_x[0]);
  a_p_abc[1] += tmp*a_y[2];
  tmp = (a_x[3]-a_x[2]) * (a_x[3]-a_x[0]) * (a_x[2]-a_x[0]) * (a_x[3]+a_x[2]-2*a_x[0]);
  a_p_abc[1] -= tmp*a_y[1];
  tmp = (a_x[3]-a_x[2]) * (a_x[3]-a_x[1]) * (a_x[2]-a_x[1]) * (a_x[3]+a_x[2]+a_x[1]-3*a_x[0]);
  a_p_abc[1] += tmp*a_y[0];
  
  /// Doing "c" term by term
  tmp = (a_x[2]-a_x[1]) * (a_x[2]-a_x[0]) * (a_x[1]-a_x[0]);
  a_p_abc[2] = tmp*a_y[3];
  tmp = (a_x[3]-a_x[1]) * (a_x[3]-a_x[0]) * (a_x[1]-a_x[0]);
  a_p_abc[2] -= tmp*a_y[2];
  tmp = (a_x[3]-a_x[2]) * (a_x[3]-a_x[0]) * (a_x[2]-a_x[0]);
  a_p_abc[2] += tmp*a_y[1];
  tmp = (a_x[3]-a_x[2]) * (a_x[3]-a_x[1]) * (a_x[2]-a_x[1]);
  a_p_abc[2] -= tmp*a_y[0];
  
  /// Normalization to D
  a_p_abc[0] = a_p_abc[0]/D3/D2;
  a_p_abc[1] = a_p_abc[1]/D3/D2;
  a_p_abc[2] = a_p_abc[2]/D3/D2;
  
  /// Prepare output
  tmp = a_y[0] + (a_x_f-a_x[0]) * a_p_abc[0];
  tmp += (a_x_f-a_x[0]) * (a_x_f-a_x[0]) * a_p_abc[1];
  tmp += (a_x_f-a_x[0]) * (a_x_f-a_x[0]) * (a_x_f-a_x[0])*a_p_abc[2];
  /// Change the polynomial to x_f-centered: p_abc remains the same, p_abc[1] gets an extra 3*c*(x_f-x[0])
  a_p_abc[1] += 3*a_p_abc[2]*(a_x_f-a_x[0]);
  /// p_abc[0] written in the new p_abc[1], so it has minus sign
  a_p_abc[0] += 2*a_p_abc[1] * (a_x_f-a_x[0]);
  a_p_abc[0] -= 3*a_p_abc[2] * (a_x_f-a_x[0]) * (a_x_f-a_x[0]);
  return tmp;
}

bool threeDiagSolver(double* a_x, double* a_a, double* a_b, double* a_c, double* a_f, int a_N, double a_ERR)
{
  /** A standard method of solving 3-diagonal matrix equation in N steps.
   *  Assume we have A*X=F, where A is NxN matrix with three diagonals.
   *  The main diagonal is a_0 .. a_{N-1}, the right diagonal b_0 .. B_{N-2}, left diagonal is c_0 .. c_{N-2},
   *  so that c_{i-1}*x_{i-1} + a_i*x_i + b_i*x_{x+1} = f_i.
   *  The assumption is that x_i = \alpha_i*x_{i+1} + \beta_i,
   *  where \alpha_i and \beta_i are unknowns.
   *  From the first line we obtain
   *  \alpha_0 = -b_0 / a_0;
   *  \beta_0 = f_0 / a_0;
   *  and for all others
   *  \alpha_i = -b_i / (a_i + c_{i-1}*\alpha_{i-1});
   *  \beta_i = (f_i - c_{i-1}*\beta_{i-1}) / (a_i + c_{i-1}*\alpha_{i-1}).
   *  \alpha and \beta go from _0 to _{N-2}.
   *  As the result, vector x is the solution and all other vectors and matrices are CHANGED!
   *  Function returns TRUE if successful, FALSE otherwise.
   *  The solver does not always solve the system even when Det(M)!=0.
   *  However, there is matrix reshuffling if the solver breaks at the first step.
   */

  if(a_N<2) {std::cerr << "ERROR: threeDiagSolver() failed! Too small matrix size!\n"; return false;}
  bool swap_lines = false;
  if (fabs(a_a[0])<a_ERR){
    if(fabs(a_b[0])<a_ERR) {std::cerr << "ERROR: threeDiagSolver() failed! Determinant is zero.\n"; return false;}
    /// Attempt to reshuffle matrix
    swap_lines = true;
    a_a[0] = a_b[0];
    double y0 = a_f[0]/a_a[0];
    a_b[0] = 0;
    a_f[1] = a_f[1] - y0*a_a[1];
    a_a[1] = a_c[0];
    a_c[0] = 0;
    a_f[2] = a_f[2] - y0*a_c[1];
    a_c[1] = 0;
  }
  double* alpha = new double[a_N-1];
  double* beta = new double[a_N-1];
  
  /// Create alpha_i and beta_i
  alpha[0] = -a_b[0]/a_a[0];
  beta[0] = a_f[0]/a_a[0];
  double D;
  
  for(int i=1;i<a_N-1;i++){
    D = (a_a[i] + a_c[i-1]*alpha[i-1]);
    if (fabs(D)<a_ERR) {std::cerr << "ERROR: threeDiagSolver() failed! a[" << i << "] is too small.\n"; return false;}
    alpha[i] = -a_b[i]/D;
    beta[i] = (a_f[i] - a_c[i-1]*beta[i-1])/D;
  }
  
  /// Sweep backward
  D = (a_a[a_N-1] + a_c[a_N-2]*alpha[a_N-2]);
  if (fabs(D)<a_ERR) {std::cerr << "ERROR: threeDiagSolver() failed at the backward sweep stage!\n"; return false;}
  a_x[a_N-1] = (a_f[a_N-1] - a_c[a_N-2]*beta[a_N-2])/D;
  
  for(int i=2;i<a_N+1;i++){
    a_x[a_N-i] = alpha[a_N-i]*a_x[a_N-i+1] + beta[a_N-i];
  }
  if(swap_lines){
    double tmp = a_x[0];
    a_x[0] = a_x[1];
    a_x[1] = tmp;
  }
  
  delete[] alpha;
  delete[] beta;
  return true;
}

#include "NamespaceFooter.H"
