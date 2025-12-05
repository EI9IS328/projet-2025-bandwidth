//************************************************************************
//   proxy application v.0.0.1
//
//  semproxy.cpp: the main interface of  proxy application
//
//************************************************************************

#include "sem_proxy.h"

#include <cartesian_struct_builder.h>
#include <cartesian_unstruct_builder.h>
#include <sem_solver_acoustic.h>
#include <source_and_receiver_utils.h>

#include <cxxopts.hpp>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <variant>

using namespace SourceAndReceiverUtils;

SEMproxy::SEMproxy(const SemProxyOptions& opt)
{
  const int order = opt.order;
  nb_elements_[0] = opt.ex;
  nb_elements_[1] = opt.ey;
  nb_elements_[2] = opt.ez;
  nb_nodes_[0] = opt.ex * order + 1;
  nb_nodes_[1] = opt.ey * order + 1;
  nb_nodes_[2] = opt.ez * order + 1;
  snapshot = opt.snapshot; // ADDED the SNAPSHOT HERE 
  const float spongex = opt.boundaries_size;
  const float spongey = opt.boundaries_size;
  const float spongez = opt.boundaries_size;
  const std::array<float, 3> sponge_size = {spongex, spongey, spongez};
  src_coord_[0] = opt.srcx;
  src_coord_[1] = opt.srcy;
  src_coord_[2] = opt.srcz;

  domain_size_[0] = opt.lx;
  domain_size_[1] = opt.ly;
  domain_size_[2] = opt.lz;

  rcv_coord_[0] = opt.rcvx;
  rcv_coord_[1] = opt.rcvy;
  rcv_coord_[2] = opt.rcvz;

  bool isModelOnNodes = opt.isModelOnNodes;
  isElastic_ = opt.isElastic;
  cout << boolalpha;
  bool isElastic = isElastic_;

  const SolverFactory::methodType methodType = getMethod(opt.method);
  const SolverFactory::implemType implemType = getImplem(opt.implem);
  const SolverFactory::meshType meshType = getMesh(opt.mesh);
  const SolverFactory::modelLocationType modelLocation =
      isModelOnNodes ? SolverFactory::modelLocationType::OnNodes
                     : SolverFactory::modelLocationType::OnElements;
  const SolverFactory::physicType physicType = SolverFactory::physicType::Acoustic;

  float lx = domain_size_[0];
  float ly = domain_size_[1];
  float lz = domain_size_[2];
  int ex = nb_elements_[0];
  int ey = nb_elements_[1];
  int ez = nb_elements_[2];

  if (meshType == SolverFactory::Struct)
  {
    switch (order)
    {
      case 1: {
        model::CartesianStructBuilder<float, int, 1> builder(
            ex, lx, ey, ly, ez, lz, isModelOnNodes);
        m_mesh = builder.getModel();
        break;
      }
      case 2: {
        model::CartesianStructBuilder<float, int, 2> builder(
            ex, lx, ey, ly, ez, lz, isModelOnNodes);
        m_mesh = builder.getModel();
        break;
      }
      case 3: {
        model::CartesianStructBuilder<float, int, 3> builder(
            ex, lx, ey, ly, ez, lz, isModelOnNodes);
        m_mesh = builder.getModel();
        break;
      }
      default:
        throw std::runtime_error(
            "Order other than 1 2 3 is not supported (semproxy)");
    }
  }
  else if (meshType == SolverFactory::Unstruct)
  {
    model::CartesianParams<float, int> param(order, ex, ey, ez, lx, ly, lz,
                                             isModelOnNodes);
    model::CartesianUnstructBuilder<float, int> builder(param);
    m_mesh = builder.getModel();
  }
  else
  {
    throw std::runtime_error("Incorrect mesh type (SEMproxy ctor.)");
  }

  // time parameters
  if (opt.autodt)
  {
    float cfl_factor = (order == 2) ? 0.5 : 0.7;
    dt_ = find_cfl_dt(cfl_factor);
  }
  else
  {
    dt_ = opt.dt;
  }
  timemax_ = opt.timemax;
  num_sample_ = timemax_ / dt_;

  m_solver = SolverFactory::createSolver(methodType, implemType, meshType,
                                         modelLocation, physicType, order);
  m_solver->computeFEInit(*m_mesh, sponge_size, opt.surface_sponge,
                          opt.taper_delta);

  initFiniteElem();

  std::cout << "Number of node is " << m_mesh->getNumberOfNodes() << std::endl;
  std::cout << "Number of element is " << m_mesh->getNumberOfElements()
            << std::endl;
  std::cout << "Launching the Method " << opt.method << ", the implementation "
            << opt.implem << " and the mesh is " << opt.mesh << std::endl;
  std::cout << "Model is on " << (isModelOnNodes ? "nodes" : "elements")
            << std::endl;
  std::cout << "Physics type is " << (isElastic ? "elastic" : "acoustic")
            << std::endl;
  std::cout << "Order of approximation will be " << order << std::endl;
  std::cout << "Time step is " << dt_ << "s" << std::endl;
  std::cout << "Simulated time is " << timemax_ << "s" << std::endl;

}

void SEMproxy::run()
{
  time_point<system_clock> startComputeTime, startOutputTime, totalComputeTime,
      totalOutputTime;

  SEMsolverDataAcoustic solverData(i1, i2, myRHSTerm, pnGlobal, rhsElement,
                                   rhsWeights);
  	std::ofstream out;
    if (snapshot > 0)
        out.open("results.txt", std::ios::app);
    

    if (snapshot > 0 && !out)
        std::cerr << "Error: cannot open results.txt\n";
    // --- 🆕 AJOUT : écrire les paramètres au début du fichier ---
    if (snapshot > 0 && out) {
        out << "# Simulation parameters\n";
        out << "# ex=" << nb_elements_[0] << "\n";
        out << "# ey=" << nb_elements_[1] << "\n";
        out << "# ez=" << nb_elements_[2] << "\n";
        out << "# order=" << m_mesh->getOrder() << "\n";
        out << "# dt=" << dt_ << "\n";
        out << "# timemax=" << timemax_ << "\n";
        out << "# src=(" << src_coord_[0] << "," << src_coord_[1] << "," << src_coord_[2] << ")\n";
        out << "# rcv=(" << rcv_coord_[0] << "," << rcv_coord_[1] << "," << rcv_coord_[2] << ")\n";
        out << "# Columns: Step, Elem, i, j, k, X, Y, Z, pnGlobal\n\n";
    }
  for (int indexTimeSample = 0; indexTimeSample < num_sample_;
       indexTimeSample++)
  {
    startComputeTime = system_clock::now();
    m_solver->computeOneStep(dt_, indexTimeSample, solverData);
    totalComputeTime += system_clock::now() - startComputeTime;

    startOutputTime = system_clock::now();

    if (indexTimeSample % 50 == 0)
    {
      m_solver->outputSolutionValues(indexTimeSample, i1, rhsElement[0],
                                     pnGlobal, "pnGlobal");
    }
    // --- SNAPSHOT OUTPUT: Step x y z pnGlobal ---
    if (snapshot > 0 &&
        indexTimeSample != 0 &&
        indexTimeSample % snapshot == 0 &&
        out)  // fichier bien ouvert
    {
      int ne    = m_mesh->getNumberOfElements();
      int order = m_mesh->getOrder();

      for (int e = 0; e < ne; ++e) {
        for (int k = 0; k <= order; ++k) {
          for (int j = 0; j <= order; ++j) {
            for (int i = 0; i <= order; ++i) {

              int nodeIdx = m_mesh->globalNodeIndex(e, i, j, k);

              float x = m_mesh->nodeCoord(nodeIdx, 0);
              float y = m_mesh->nodeCoord(nodeIdx, 1);
              float z = m_mesh->nodeCoord(nodeIdx, 2);

              float value = pnGlobal(nodeIdx, i1);

              // Step  elem  i  j  k  x  y  z  pnGlobal
              out << indexTimeSample << " "
                  << e   << " "
                  << i   << " " << j << " " << k << " "
                  << x   << " " << y << " " << z << " "
                  << value << "\n";
            }
          }
        }
      }
    }

    
// Save pressure at receiver
    float varnp1 = 0.0;
    const int order = m_mesh->getOrder();
    for(int m = 0;m< num_receivers;m++){
      varnp1 = pnGlobal(rhsElementRcv[m], i2);
      float xn = m_mesh->nodeCoord(rhsElementRcv[m], 0);
      float yn = m_mesh->nodeCoord(rhsElementRcv[m], 1);
      float zn = m_mesh->nodeCoord(rhsElementRcv[m], 2);

      std::ofstream recv("recev_results_" + std::to_string(m) + ".txt", std::ios::app);
      recv << indexTimeSample << " " << xn << " " << yn << " " << zn << " " << varnp1 << "\n";
      pnAtReceiver(m, indexTimeSample) = varnp1;
      }

    swap(i1, i2);

    auto tmp = solverData.m_i1;
    solverData.m_i1 = solverData.m_i2;
    solverData.m_i2 = tmp;

    totalOutputTime += system_clock::now() - startOutputTime;
  }

  float kerneltime_ms = time_point_cast<microseconds>(totalComputeTime)
                            .time_since_epoch()
                            .count();
  float outputtime_ms =
      time_point_cast<microseconds>(totalOutputTime).time_since_epoch().count();

  cout << "------------------------------------------------ " << endl;
  cout << "\n---- Elapsed Kernel Time : " << kerneltime_ms / 1E6 << " seconds."
       << endl;
  cout << "---- Elapsed Output Time : " << outputtime_ms / 1E6 << " seconds."
       << endl;
  cout << "------------------------------------------------ " << endl;
}

bool pointInDomain(float x, float y, float z,
                   float width, float height, float length)
{
    return (x >= 0 && x <= width &&
            y >= 0 && y <= height &&
            z >= 0 && z <= length);
}

int countValidPoints(const std::string& filename,
                     float width, float height, float length)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Erreur : impossible d'ouvrir le fichier : " << filename << std::endl;
        return 1;
    }

    float x, y, z;
    int validCount = 0;
    int lineNumber = 0;

    while (file >> x >> y >> z) {   // Lit 3 valeurs par ligne
        lineNumber++;

        if (pointInDomain(x, y, z, width, height, length)) {
            validCount++;
        }
    }

    return validCount;
}


// Initialize arrays
void SEMproxy::init_arrays()
{
  cout << "Allocate host memory for source and pressure values ..." << endl;
  num_receivers = countValidPoints(recv_file,domain_size_[0],domain_size_[1],domain_size_[2]);
  rhsElement = allocateVector<vectorInt>(myNumberOfRHS, "rhsElement");
  rhsWeights = allocateArray2D<arrayReal>(
      myNumberOfRHS, m_mesh->getNumberOfPointsPerElement(), "RHSWeight");
  myRHSTerm = allocateArray2D<arrayReal>(myNumberOfRHS, num_sample_, "RHSTerm");
  pnGlobal =
      allocateArray2D<arrayReal>(m_mesh->getNumberOfNodes(), 2, "pnGlobal");
  pnAtReceiver = allocateArray2D<arrayReal>(num_receivers, num_sample_, "pnAtReceiver");
  // Receiver
  rhsElementRcv = allocateVector<vectorInt>(num_receivers, "rhsElementRcv");
  rhsWeightsRcv = allocateArray2D<arrayReal>(
      num_receivers, m_mesh->getNumberOfPointsPerElement(), "RHSWeightRcv");
}


float intra_distance(float xn ,float yn ,float zn, float x, float y, float z)
{
    return std::sqrt(
        (x - xn) * (x - xn) +
        (y - yn) * (y - yn) +
        (z - zn) * (z - zn)
    );
}


// Initialize sources
void SEMproxy::init_source()
{
  arrayReal myRHSLocation = allocateArray2D<arrayReal>(1, 3, "RHSLocation");
  // std::cout << "All source are currently are coded on element 50." <<
  // std::endl;
  std::cout << "All source are currently are coded on middle element."
            << std::endl;
  int ex = nb_elements_[0];
  int ey = nb_elements_[1];
  int ez = nb_elements_[2];

  int lx = domain_size_[0];
  int ly = domain_size_[1];
  int lz = domain_size_[2];

  // Get source element index

  int source_index = floor((src_coord_[0] * ex) / lx) +
                     floor((src_coord_[1] * ey) / ly) * ex +
                     floor((src_coord_[2] * ez) / lz) * ey * ex;

  for (int i = 0; i < 1; i++)
  {
    rhsElement[i] = source_index;
  }

  // Get coordinates of the corners of the sourc element
  float cornerCoords[8][3];
  int I = 0;
  int nodes_corner[2] = {0, m_mesh->getOrder()};
  for (int k : nodes_corner)
  {
    for (int j : nodes_corner)
    {
      for (int i : nodes_corner)
      {
        int nodeIdx = m_mesh->globalNodeIndex(rhsElement[0], i, j, k);
        cornerCoords[I][0] = m_mesh->nodeCoord(nodeIdx, 0);
        cornerCoords[I][2] = m_mesh->nodeCoord(nodeIdx, 2);
        cornerCoords[I][1] = m_mesh->nodeCoord(nodeIdx, 1);
        I++;
      }
    }
  }

  // initialize source term
  vector<float> sourceTerm =
      myUtils.computeSourceTerm(num_sample_, dt_, f0, sourceOrder);
  for (int j = 0; j < num_sample_; j++)
  {
    myRHSTerm(0, j) = sourceTerm[j];
    if (j % 100 == 0)
      cout << "Sample " << j << "\t: sourceTerm = " << sourceTerm[j] << endl;
  }

  // get element number of source term
  myElementSource = rhsElement[0];
  cout << "Element number for the source location: " << myElementSource << endl
       << endl;

  int order = m_mesh->getOrder();

  switch (order)
  {
    case 1:
      SourceAndReceiverUtils::ComputeRHSWeights<1>(cornerCoords, src_coord_,
                                                   rhsWeights);
      break;
    case 2:
      SourceAndReceiverUtils::ComputeRHSWeights<2>(cornerCoords, src_coord_,
                                                   rhsWeights);
      break;
    case 3:
      SourceAndReceiverUtils::ComputeRHSWeights<3>(cornerCoords, src_coord_,
                                                   rhsWeights);
      break;
    default:
      throw std::runtime_error("Unsupported order: " + std::to_string(order));
  }

  // Receiver computation
  int receiver_index = floor((rcv_coord_[0] * ex) / lx) +
                       floor((rcv_coord_[1] * ey) / ly) * ex +
                       floor((rcv_coord_[2] * ez) / lz) * ey * ex;

  for (int i = 0; i < 1; i++)
  {
    rhsElementRcv[i] = receiver_index;
  }

  // Get coordinates of the corners of the receiver element
  float cornerCoordsRcv[8][3];
  I = 0;
  for (int k : nodes_corner)
  {
    for (int j : nodes_corner)
    {
      for (int i : nodes_corner)
      {
        int nodeIdx = m_mesh->globalNodeIndex(rhsElementRcv[0], i, j, k);
        cornerCoordsRcv[I][0] = m_mesh->nodeCoord(nodeIdx, 0);
        cornerCoordsRcv[I][2] = m_mesh->nodeCoord(nodeIdx, 2);
        cornerCoordsRcv[I][1] = m_mesh->nodeCoord(nodeIdx, 1);
        I++;
      }
    }
  }

//  switch (order)
//  {
//    case 1:
//      SourceAndReceiverUtils::ComputeRHSWeights<1>(cornerCoordsRcv, rcv_coord_,
//                                                   rhsWeightsRcv);
//      break;
//    case 2:
//      SourceAndReceiverUtils::ComputeRHSWeights<2>(cornerCoordsRcv, rcv_coord_,
//                                                   rhsWeightsRcv);
//      break;
//    case 3:
//      SourceAndReceiverUtils::ComputeRHSWeights<3>(cornerCoordsRcv, rcv_coord_,
//                                                   rhsWeightsRcv);
//      break;
//    default:
//      throw std::runtime_error("Unsupported order: " + std::to_string(order));
//  }

  std::ifstream file(recv_file);
    if (!file.is_open()) {
        std::cerr << "Erreur : impossible d'ouvrir le fichier : " << recv_file << std::endl;
        // Receiver computation
        int receiver_index = floor((rcv_coord_[0] * ex) / lx) +
                       floor((rcv_coord_[1] * ey) / ly) * ex +
                       floor((rcv_coord_[2] * ez) / lz) * ey * ex;
        
        rhsElementRcv[0]= receiver_index;
    }

    else {
      float x, y, z;
      int validCount = 0;
      int lineNumber = 0;
      float width = domain_size_[0],height = domain_size_[1],length = domain_size_[2];

      while (file >> x >> y >> z) {   // Lit 3 valeurs par ligne
          lineNumber++;
          int receiver_index = floor((rcv_coord_[0] * ex) / lx) +
                       floor((rcv_coord_[1] * ey) / ly) * ex +
                       floor((rcv_coord_[2] * ez) / lz) * ey * ex;
        
          rhsElementRcv[validCount]= receiver_index;

          if (pointInDomain(x, y, z, width, height, length)) {
              float cornerCoordsRcv[3];
              float distance = std::numeric_limits<float>::infinity();
              for (int k : nodes_corner)
              {
                for (int j : nodes_corner)
                {
                  for (int i : nodes_corner)
                  {
                    int nodeIdx = m_mesh->globalNodeIndex(rhsElementRcv[validCount], i, j, k);
                    float xn = m_mesh->nodeCoord(nodeIdx, 0);
                    float yn = m_mesh->nodeCoord(nodeIdx, 1);
                    float zn = m_mesh->nodeCoord(nodeIdx, 2);
                    if(intra_distance(xn,yn,zn,x,y,z) < distance){
                      receiver_index = nodeIdx;
                    }
                  }
                }
              }

              rhsElementRcv[validCount]= receiver_index;
              validCount++;
          }
      }
    }
}

SolverFactory::implemType SEMproxy::getImplem(string implemArg)
{
  if (implemArg == "makutu") return SolverFactory::MAKUTU;
  if (implemArg == "shiva") return SolverFactory::SHIVA;

  throw std::invalid_argument(
      "Implentation type does not follow any valid type.");
}

SolverFactory::meshType SEMproxy::getMesh(string meshArg)
{
  if (meshArg == "cartesian") return SolverFactory::Struct;
  if (meshArg == "ucartesian") return SolverFactory::Unstruct;

  std::cout << "Mesh type found is " << meshArg << std::endl;
  throw std::invalid_argument("Mesh type does not follow any valid type.");
}

SolverFactory::methodType SEMproxy::getMethod(string methodArg)
{
  if (methodArg == "sem") return SolverFactory::SEM;
  if (methodArg == "dg") return SolverFactory::DG;

  throw std::invalid_argument("Method type does not follow any valid type.");
}

float SEMproxy::find_cfl_dt(float cfl_factor)
{
  float sqrtDim3 = 1.73;  // to change for 2d
  float min_spacing = m_mesh->getMinSpacing();
  float v_max = m_mesh->getMaxSpeed();

  float dt = cfl_factor * min_spacing / (sqrtDim3 * v_max);

  return dt;
}
