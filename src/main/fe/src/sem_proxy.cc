//************************************************************************
//   proxy application v.0.0.1
//
//  semproxy.cpp: the main interface of  proxy application
//
//************************************************************************

#include <cartesian_struct_builder.h>
#include <cartesian_unstruct_builder.h>
#include <model.h>
#include <sem_solver_acoustic.h>
#include <source_and_receiver_utils.h>

#include <cxxopts.hpp>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <variant>

#include "sem_proxy.h"

using namespace SourceAndReceiverUtils;

struct RGB
{
  unsigned char r, g, b;
};

RGB colormap(float v)
{
  v = std::max(-1.0f, std::min(1.0f, v));
  if (v < 0)
    return {0, static_cast<unsigned char>(255 * (1 + v)), 255};
  else
    return {255, static_cast<unsigned char>(255 * (1 - v)), 0};
}

void writePPMSliceXY_single(const std::string& filename, arrayReal& pnGlobal,
                            model::ModelApi<float, int>* mesh)
{
  int order = mesh->getOrder();
  int nx = mesh->getNumberOfElements() * order + 1;
  int ny = mesh->getNumberOfElements() * order + 1;

  int kSlice = mesh->getNumberOfElements() * order / 2;  // slice centrale

  std::ofstream ppm(filename, std::ios::binary);
  ppm << "P6\n" << nx << " " << ny << "\n255\n";

  float minV = 1e30f, maxV = -1e30f;

  // trouver min et max sur la slice
  for (int j = 0; j < ny; j++)
    for (int i = 0; i < nx; i++)
    {
      int node = mesh->globalNodeIndex(0, i, j, kSlice);
      float v = pnGlobal(node, 0);  // ici timestep 0 ou dernier timestep
      minV = std::min(minV, v);
      maxV = std::max(maxV, v);
    }

  float invRange = (maxV > minV) ? 1.0f / (maxV - minV) : 1.0f;

  // écrire la slice
  for (int j = 0; j < ny; j++)
  {
    for (int i = 0; i < nx; i++)
    {
      int node = mesh->globalNodeIndex(0, i, j, kSlice);
      float v = pnGlobal(node, 0);
      float vn = 2.0f * (v - minV) * invRange - 1.0f;
      RGB c = colormap(vn);
      ppm.write(reinterpret_cast<char*>(&c), 3);
    }
  }
}

// Ajout de dequantizeSnapshot
void dequantizeSnapshot(const std::vector<uint16_t>& input,
                        std::vector<float>& output, float minV, float maxV,
                        int nbits)
{
  int levels = (1 << nbits) - 1;
  output.resize(input.size());
  for (size_t i = 0; i < input.size(); i++)
  {
    float norm = static_cast<float>(input[i]) / levels;
    output[i] = minV + norm * (maxV - minV);
  }
}

// ====================== COMPRESSION ===========================

float computeRMSE(const std::vector<float>& a, const std::vector<float>& b)
{
  float err = 0.0f;
  for (size_t i = 0; i < a.size(); i++) err += (a[i] - b[i]) * (a[i] - b[i]);
  return std::sqrt(err / a.size());
}

void quantizeSnapshot(const std::vector<float>& input,
                      std::vector<uint16_t>& output, float& minV, float& maxV,
                      int nbits)
{
  minV = *std::min_element(input.begin(), input.end());
  maxV = *std::max_element(input.begin(), input.end());

  int levels = (1 << nbits) - 1;
  output.resize(input.size());

  for (size_t i = 0; i < input.size(); i++)
  {
    float norm = (input[i] - minV) / (maxV - minV);
    output[i] = static_cast<uint16_t>(norm * levels);
  }
}

std::vector<std::pair<uint16_t, uint16_t>> RLEencode(
    const std::vector<uint16_t>& data)
{
  std::vector<std::pair<uint16_t, uint16_t>> out;

  uint16_t val = data[0], count = 1;
  for (size_t i = 1; i < data.size(); i++)
  {
    if (data[i] == val && count < 65535)
      count++;
    else
    {
      out.emplace_back(val, count);
      val = data[i];
      count = 1;
    }
  }
  out.emplace_back(val, count);
  return out;
}

SEMproxy::SEMproxy(const SemProxyOptions& opt)
{
  const int order = opt.order;
  nb_elements_[0] = opt.ex;
  nb_elements_[1] = opt.ey;
  nb_elements_[2] = opt.ez;
  nb_nodes_[0] = opt.ex * order + 1;
  nb_nodes_[1] = opt.ey * order + 1;
  nb_nodes_[2] = opt.ez * order + 1;
  snapshot = opt.snapshot;  // ADDED the SNAPSHOT HERE
  recv_file = opt.recv_file;
  recv_on = opt.recv_on;
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
  const SolverFactory::physicType physicType =
      SolverFactory::physicType::Acoustic;

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
  using std::chrono::duration_cast;
  using std::chrono::microseconds;
  using std::chrono::system_clock;
  using std::chrono::time_point;

  // Timers for compute and output
  time_point<system_clock> startComputeTime, startOutputTime;
  system_clock::duration totalComputeTime = system_clock::duration::zero();
  system_clock::duration totalOutputTime = system_clock::duration::zero();

  // Extra accumulators (in microseconds) for detailed breakdown
  double snapshotTime_us = 0.0;
  double sismoTime_us = 0.0;

  // Measure total wall-clock of run()
  auto totalStart = system_clock::now();

  SEMsolverDataAcoustic solverData(i1, i2, myRHSTerm, pnGlobal, rhsElement,
                                   rhsWeights);

  // Snapshot file (volumetric output)
  std::ofstream out;
  if (snapshot > 0)
  {
    out.open("results.csv", std::ios::app);
    if (!out)
    {
      std::cerr << "Error: cannot open results.csv\n";
    }
    else
    {
      if (out.tellp() == 0)
      {
        out << "Step,X,Y,Z,pnGlobal\n";
      }
    }
  }

  for (int indexTimeSample = 0; indexTimeSample < num_sample_;
       indexTimeSample++)
  {
    // =======================
    // 1) Kernel / compute
    // =======================
    startComputeTime = system_clock::now();
    m_solver->computeOneStep(dt_, indexTimeSample, solverData);
    totalComputeTime += system_clock::now() - startComputeTime;

    // =======================
    // 2) Output region
    // =======================
    startOutputTime = system_clock::now();

    if (indexTimeSample % 50 == 0)
    {
      m_solver->outputSolutionValues(indexTimeSample, i1, rhsElement[0],
                                     pnGlobal, "pnGlobal");
    }

    // Après la boucle for (int indexTimeSample = 0; ...) // run()
    if (snapshot > 0 && )
    {
      writePPMSliceXY_single("slice_xy.ppm", pnGlobal, m_mesh.get());
    }

    // --- SNAPSHOT timing + writing ---
    if (snapshot > 0 && indexTimeSample != 0 &&
        indexTimeSample % snapshot == 0 && out)  // file OK
    {
      auto startSnapshot = system_clock::now();

      int ne = m_mesh->getNumberOfElements();
      int order = m_mesh->getOrder();

      for (int e = 0; e < ne; ++e)
      {
        for (int k = 0; k <= order; ++k)
        {
          for (int j = 0; j <= order; ++j)
          {
            for (int i = 0; i <= order; ++i)
            {
              int nodeIdx = m_mesh->globalNodeIndex(e, i, j, k);

              float x = m_mesh->nodeCoord(nodeIdx, 0);
              float y = m_mesh->nodeCoord(nodeIdx, 1);
              float z = m_mesh->nodeCoord(nodeIdx, 2);

              float value = pnGlobal(nodeIdx, i1);

              out << indexTimeSample << "," << x << "," << y << "," << z << ","
                  << value << "\n";
            }
          }
        }
      }

      auto endSnapshot = system_clock::now();
      snapshotTime_us +=
          duration_cast<microseconds>(endSnapshot - startSnapshot).count();
    }

    // ===== TP3 : Quantification + RMSE =====
    if (snapshot > 0 && indexTimeSample != 0 && indexTimeSample % snapshot == 0)
    {
      std::vector<float> original;
      std::vector<float> reconstructed;
      std::vector<uint16_t> quantized;

      int nbNodes = m_mesh->getNumberOfNodes();
      original.reserve(nbNodes);

      for (int n = 0; n < nbNodes; n++) original.push_back(pnGlobal(n, i1));

      float minV, maxV;
      int nbits = 8;  // <- test 8 / 16 / 32 bits

      quantizeSnapshot(original, quantized, minV, maxV, nbits);
      dequantizeSnapshot(quantized, reconstructed, minV, maxV, nbits);

      float rmse = computeRMSE(original, reconstructed);

      std::ofstream errFile("quantization_error.csv", std::ios::app);
      if (errFile.tellp() == 0) errFile << "timeStep,nbits,rmse\n";

      errFile << indexTimeSample << "," << nbits << "," << rmse << "\n";
    }

    // =======================
    // 3) Sismo / receivers
    // =======================
    if (recv_on && num_receivers > 0)
    {
      auto startSismo = system_clock::now();
      float varnp1 = 0.0f;

      for (int m = 0; m < num_receivers; m++)
      {
        varnp1 = pnGlobal(rhsElementRcv[m], i2);

        float xn = m_mesh->nodeCoord(rhsElementRcv[m], 0);
        float yn = m_mesh->nodeCoord(rhsElementRcv[m], 1);
        float zn = m_mesh->nodeCoord(rhsElementRcv[m], 2);

        if (indexTimeSample == 0)
        {
          std::ofstream recv("recev_results_" + std::to_string(m) + ".csv",
                             std::ios::app);
          if (recv.tellp() == 0) recv << "indexTimeSample,xn,yn,zn,varnp1\n";
        }

        std::ofstream recv("recev_results_" + std::to_string(m) + ".csv",
                           std::ios::app);
        recv << indexTimeSample << "," << xn << "," << yn << "," << zn << ","
             << varnp1 << "\n";

        pnAtReceiver(m, indexTimeSample) = varnp1;
      }

      auto endSismo = system_clock::now();
      sismoTime_us +=
          duration_cast<microseconds>(endSismo - startSismo).count();
    }

    // =======================
    // 4) Swap buffers
    // =======================
    std::swap(i1, i2);

    auto tmp = solverData.m_i1;
    solverData.m_i1 = solverData.m_i2;
    solverData.m_i2 = tmp;

    totalOutputTime += system_clock::now() - startOutputTime;
  }

  auto totalEnd = system_clock::now();
  double totalRun_us =
      duration_cast<microseconds>(totalEnd - totalStart).count();

  double kerneltime_us = duration_cast<microseconds>(totalComputeTime).count();
  double outputtime_us = duration_cast<microseconds>(totalOutputTime).count();

  double kernel_s = kerneltime_us / 1e6;
  double output_s = outputtime_us / 1e6;
  double snapshot_s = snapshotTime_us / 1e6;
  double sismo_s = sismoTime_us / 1e6;
  double total_s = totalRun_us / 1e6;

  std::cout << "------------------------------------------------ \n";
  std::cout << "---- Elapsed Kernel Time   : " << kernel_s << " s\n";
  std::cout << "---- Elapsed Output Time   : " << output_s << " s\n";
  std::cout << "----   - Snapshot Time     : " << snapshot_s << " s\n";
  std::cout << "----   - Sismo Time        : " << sismo_s << " s\n";
  std::cout << "---- Total run() Time      : " << total_s << " s\n";
  std::cout << "------------------------------------------------ \n";

  // SAVE RUN TIME INFORMATION TO CSV
  std::ofstream tfile("time_debug.csv", std::ios::app);
  if (tfile)
  {
    if (tfile.tellp() == 0)
      tfile
          << "ex,ey,ez,snapshot,kernel_s,output_s,snapshot_s,sismo_s,total_s\n";

    int ex = nb_elements_[0];
    int ey = nb_elements_[1];
    int ez = nb_elements_[2];

    tfile << ex << "," << ey << "," << ez << "," << snapshot << "," << kernel_s
          << "," << output_s << "," << snapshot_s << "," << sismo_s << ","
          << total_s << "\n";
  }
}

bool pointInDomain(float x, float y, float z, float width, float height,
                   float length)
{
  return (x >= 0 && x <= width && y >= 0 && y <= height && z >= 0 &&
          z <= length);
}

int countValidPoints(const std::string& filename, float width, float height,
                     float length)
{
  std::ifstream file(filename);
  if (!file.is_open())
  {
    return 1;
  }

  float x, y, z;
  int validCount = 0;
  int lineNumber = 0;

  while (file >> x >> y >> z)
  {  // Lit 3 valeurs par ligne
    lineNumber++;

    if (pointInDomain(x, y, z, width, height, length))
    {
      validCount++;
    }
  }

  return validCount;
}

// Initialize arrays
void SEMproxy::init_arrays()
{
  cout << "Allocate host memory for source and pressure values ..." << recv_file
       << endl;
  num_receivers = countValidPoints(recv_file, domain_size_[0], domain_size_[1],
                                   domain_size_[2]);
  rhsElement = allocateVector<vectorInt>(myNumberOfRHS, "rhsElement");
  rhsWeights = allocateArray2D<arrayReal>(
      myNumberOfRHS, m_mesh->getNumberOfPointsPerElement(), "RHSWeight");
  myRHSTerm = allocateArray2D<arrayReal>(myNumberOfRHS, num_sample_, "RHSTerm");
  pnGlobal =
      allocateArray2D<arrayReal>(m_mesh->getNumberOfNodes(), 2, "pnGlobal");
  pnAtReceiver =
      allocateArray2D<arrayReal>(num_receivers, num_sample_, "pnAtReceiver");
  // Receiver
  rhsElementRcv = allocateVector<vectorInt>(num_receivers, "rhsElementRcv");
  rhsWeightsRcv = allocateArray2D<arrayReal>(
      num_receivers, m_mesh->getNumberOfPointsPerElement(), "RHSWeightRcv");
}

float intra_distance(float xn, float yn, float zn, float x, float y, float z)
{
  return std::sqrt((x - xn) * (x - xn) + (y - yn) * (y - yn) +
                   (z - zn) * (z - zn));
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
  cout << rcv_coord_[0] << " " << rcv_coord_[1] << endl;
  // Receiver computation
  int receiver_index = floor((rcv_coord_[0] * ex) / lx) +
                       floor((rcv_coord_[1] * ey) / ly) * ex +
                       floor((rcv_coord_[2] * ez) / lz) * ey * ex;

  float bestDist = std::numeric_limits<float>::infinity();
  int bestNode = -1;

  for (int k : nodes_corner)
    for (int j : nodes_corner)
      for (int i : nodes_corner)
      {
        int nodeIdx = m_mesh->globalNodeIndex(receiver_index, i, j, k);

        float xn = m_mesh->nodeCoord(nodeIdx, 0);
        float yn = m_mesh->nodeCoord(nodeIdx, 1);
        float zn = m_mesh->nodeCoord(nodeIdx, 2);

        float d = intra_distance(xn, yn, zn, rcv_coord_[0], rcv_coord_[1],
                                 rcv_coord_[2]);
        if (d < bestDist)
        {
          bestDist = d;
          bestNode = nodeIdx;
        }
      }

  rhsElementRcv[0] = bestNode;

  ////  switch (order)
  ////  {
  //    case 1:
  //      SourceAndReceiverUtils::ComputeRHSWeights<1>(cornerCoordsRcv,
  //      rcv_coord_,
  //                                                   rhsWeightsRcv);
  //      break;
  //    case 2:
  //      SourceAndReceiverUtils::ComputeRHSWeights<2>(cornerCoordsRcv,
  //      rcv_coord_,
  //                                                   rhsWeightsRcv);
  //      break;
  //    case 3:
  //      SourceAndReceiverUtils::ComputeRHSWeights<3>(cornerCoordsRcv,
  //      rcv_coord_,
  //                                                   rhsWeightsRcv);
  //      break;
  //    default:
  //      throw std::runtime_error("Unsupported order: " +
  //      std::to_string(order));
  //  }

  std::ifstream file(recv_file);
  if (!file.is_open()) return;

  float x, y, z;
  int validCount = 0;
  float width = domain_size_[0], height = domain_size_[1],
        length = domain_size_[2];

  while (file >> x >> y >> z)
  {
    if (!pointInDomain(x, y, z, width, height, length)) continue;

    // Element index from coordinates
    int element = floor((x * ex) / lx) + floor((y * ey) / ly) * ex +
                  floor((z * ez) / lz) * ey * ex;

    float bestDist = std::numeric_limits<float>::infinity();
    int bestNode = -1;

    for (int k : nodes_corner)
      for (int j : nodes_corner)
        for (int i : nodes_corner)
        {
          int nodeIdx = m_mesh->globalNodeIndex(element, i, j, k);

          float xn = m_mesh->nodeCoord(nodeIdx, 0);
          float yn = m_mesh->nodeCoord(nodeIdx, 1);
          float zn = m_mesh->nodeCoord(nodeIdx, 2);

          float d = intra_distance(xn, yn, zn, x, y, z);
          if (d < bestDist)
          {
            bestDist = d;
            bestNode = nodeIdx;
          }
        }

    rhsElementRcv[validCount++] = bestNode;
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
