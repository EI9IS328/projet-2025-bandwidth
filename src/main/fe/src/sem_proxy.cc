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

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cxxopts.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <utility>
#include <variant>
#include <vector>
using namespace SourceAndReceiverUtils;

struct RLEPair {
    uint32_t value; 
    uint32_t count;
};
struct QuantizationMeta
{
  float minVal;
  float maxVal;
  int nBits;
};

// Fonction de compression RLE
std::vector<RLEPair> RLE_compress(const std::vector<uint8_t>& data, size_t elementSize)
{
    std::vector<RLEPair> compressed;
    if (data.empty() || elementSize == 0) return compressed;

    size_t nbElements = data.size() / elementSize;


    auto getValue = [&](size_t index) -> uint32_t {
        uint32_t val = 0;
        std::memcpy(&val, &data[index * elementSize], elementSize);
        return val;
    };

    uint32_t current = getValue(0);
    uint32_t count = 1;

    for (size_t i = 1; i < nbElements; ++i)
    {
        uint32_t val = getValue(i);
        if (val == current)
        {
            ++count;
        }
        else
        {
            compressed.push_back({current, count});
            current = val;
            count = 1;
        }
    }
    compressed.push_back({current, count});

    return compressed;
}


// Écriture binaire compressée
void writeRLE(const std::string& filename,
              const std::vector<RLEPair>& compressed,
              const QuantizationMeta& qmeta)
{
  std::ofstream out(filename, std::ios::binary | std::ios::app);
  out.write(reinterpret_cast<const char*>(&qmeta), sizeof(qmeta));
  out.write(reinterpret_cast<const char*>(compressed.data()),
            compressed.size() * sizeof(RLEPair));
}

uint32_t quantize(float value, const QuantizationMeta& meta)
{
    float clamped = std::min(std::max(value, meta.minVal), meta.maxVal);
    float normalized = (clamped - meta.minVal) / (meta.maxVal - meta.minVal);

    uint64_t maxInt64;
    if (meta.nBits == 32)
        maxInt64 = uint64_t(std::numeric_limits<uint32_t>::max());
    else
        maxInt64 = (1ull << meta.nBits) - 1;

    return static_cast<uint32_t>(std::round(static_cast<double>(normalized) * maxInt64));
}


inline float dequantize(uint32_t qValue, const QuantizationMeta& meta)
{
  uint64_t maxInt64;
    if (meta.nBits == 32)
        maxInt64 = uint64_t(std::numeric_limits<uint32_t>::max());
    else
        maxInt64 = (1ull << meta.nBits) - 1;
  double normalized = static_cast<double>(qValue) / static_cast<double>(maxInt64);
  return meta.minVal + normalized * (meta.maxVal - meta.minVal);
}

float computeRMSE(const std::vector<float>& original,
                  const std::vector<float>& reconstructed)
{
  if (original.size() != reconstructed.size())
  {
    throw std::runtime_error("Vectors must have the same size");
  }
  float mse = 0.0f;
  for (size_t i = 0; i < original.size(); ++i)
  {
    float diff = original[i] - reconstructed[i];
    mse += diff * diff;
  }
  mse /= static_cast<float>(original.size());
  return std::sqrt(mse);
}

void saveSlicePPM(const std::string& filename,
                  const std::vector<float>& sliceValues, int nx, int ny,
                  float minVal, float maxVal)
{
  std::ofstream ppm(filename, std::ios::binary);
  ppm << "P6\n" << nx << " " << ny << "\n255\n";

  for (int j = 0; j < ny; ++j)
  {
    for (int i = 0; i < nx; ++i)
    {
      float val = sliceValues[j * nx + i];
      // Normalisation locale ou globale
      float normalized = (val - minVal) / (maxVal - minVal);
      normalized = std::clamp(normalized, 0.0f, 1.0f);

      // Exemple de colormap simple : jet-like (bleu → rouge)
      unsigned char r = static_cast<unsigned char>(255 * normalized);
      unsigned char g = 0;
      unsigned char b = static_cast<unsigned char>(255 * (1 - normalized));

      ppm << r << g << b;
    }
  }
  ppm.close();
}

void removeFileIfExists(const std::string& filename)
{
    if (std::filesystem::exists(filename))
    {
        std::filesystem::remove(filename);
    }
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
  ad_hoc = opt.ad_hoc;
  slices = opt.slices;
  RLE = opt.RLE;
  nBit=opt.nBit;
  quantification = opt.quantification;
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
  // In-situ statistics variables
  removeFileIfExists("snapshot_quant.bin");
  removeFileIfExists("snapshot_quant_rle.bin");
  removeFileIfExists("results.csv");
  for(int m=0;m<num_receivers;m++){
    removeFileIfExists("sismo_quant_" + std::to_string(m) + ".bin");
    removeFileIfExists("sismo_quant_rle_" + std::to_string(m) + ".bin");
    removeFileIfExists("recev_results_"+ std::to_string(m) +".csv");
  }
  float minVal = std::numeric_limits<float>::max();
  float maxVal = std::numeric_limits<float>::lowest();
  float minPreRecv[num_receivers];
  float maxPreRecv[num_receivers];
  float meanRecv[num_receivers];
  for (int m = 0; m < num_receivers; m++)
  {
    minPreRecv[m] = std::numeric_limits<float>::max();
    maxPreRecv[m] = std::numeric_limits<float>::lowest();
    meanRecv[m] = 0.0f;
  }
  //Pour calculer le taux de compression
  size_t originalSize = 0;
  size_t compressedSize = 0;
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
  startOutputTime = system_clock::now();
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
      // Write header only if file is empty
      if (out.tellp() == 0)
      {
        out << "Step,X,Y,Z,pnGlobal\n";
      }
    }
  }
  totalOutputTime += system_clock::now() - startOutputTime;

  for (int indexTimeSample = 0; indexTimeSample < num_sample_;
       indexTimeSample++)
  {
    // =======================
    // SLICE PARAMETERS
    // =======================
    const float z_slice = src_coord_[2];  // plan XY à la hauteur de la source
    const float T0 = 1.0f / f0;           // période dominante
    const float periods_per_slice = 5.0f;
    const int slice_every =
        std::max(1, static_cast<int>((periods_per_slice * T0) / dt_));
    const float tol = 0.5f * m_mesh->getMinSpacing();  // tolérance pour plan
    // =======================
    // SLICE OUTPUT (XY plane)
    // =======================
    // Écrire les slices uniquement à partir du step 1
    if (indexTimeSample > 0 && indexTimeSample % slice_every == 0 &&
        slices == true)
    {
      int nx = nb_nodes_[0];  // nombre de noeuds sur X
      int ny = nb_nodes_[1];  // nombre de noeuds sur Y
      std::vector<float> sliceValues(nx * ny, 0.0f);

      std::ofstream sliceFile;
      startOutputTime = system_clock::now();
      if (indexTimeSample == slice_every)
      {
        sliceFile.open("slice_xy.csv");
        sliceFile << "Step,X,Y,pnGlobal\n";
      }
      else
      {
        sliceFile.open("slice_xy.csv", std::ios::app);
      }
      totalOutputTime += system_clock::now() - startOutputTime;
      float tol = 0.5f * m_mesh->getMinSpacing();

      int idx = 0;
      for (int j = 0; j < ny; ++j)
      {
        for (int i = 0; i < nx; ++i)
        {
          float minDist = std::numeric_limits<float>::infinity();
          int bestNode = -1;

          // Recherche du noeud le plus proche du plan z_slice
          for (int n = 0; n < m_mesh->getNumberOfNodes(); ++n)
          {
            float xn = m_mesh->nodeCoord(n, 0);
            float yn = m_mesh->nodeCoord(n, 1);
            float zn = m_mesh->nodeCoord(n, 2);

            if (std::abs(zn - z_slice) < tol)
            {
              float dx =
                  xn - m_mesh->nodeCoord(i, 0);  // si structuré, i -> index X
              float dy = yn - m_mesh->nodeCoord(j, 1);  // j -> index Y
              float dist = std::sqrt(dx * dx + dy * dy);

              if (dist < minDist)
              {
                minDist = dist;
                bestNode = n;
              }
            }
          }

          if (bestNode == -1)
            throw std::runtime_error("No node found near XY plane at z_slice");

          float value = pnGlobal(bestNode, i1);
          sliceValues[idx++] = value;

          // Écriture CSV
          float x = m_mesh->nodeCoord(bestNode, 0);
          float y = m_mesh->nodeCoord(bestNode, 1);
          startOutputTime = system_clock::now();
          sliceFile << indexTimeSample << "," << x << "," << y << "," << value
                    << "\n";
          totalOutputTime += system_clock::now() - startOutputTime;
        }
      }

      sliceFile.close();

      float sliceMin =
          *std::min_element(sliceValues.begin(), sliceValues.end());
      float sliceMax =
          *std::max_element(sliceValues.begin(), sliceValues.end());
      startOutputTime = system_clock::now();
      saveSlicePPM("slice_xy_" + std::to_string(indexTimeSample) + ".ppm",
                   sliceValues, nx, ny, sliceMin, sliceMax);
      totalOutputTime += system_clock::now() - startOutputTime;
    }
    //=======================
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
    totalOutputTime += system_clock::now() - startOutputTime;
    // --- SNAPSHOT timing + writing ---
    if (snapshot > 0 && indexTimeSample != 0 &&
        indexTimeSample % snapshot == 0 && out)  // file OK
    {
      auto startSnapshot = system_clock::now();

      if (ad_hoc)
      {
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
                startOutputTime = system_clock::now();
                // CSV: Step, X, Y, Z, pnGlobal
                out << indexTimeSample << "," << x << "," << y << "," << z
                    << "," << value << "\n";
                totalOutputTime += system_clock::now() - startOutputTime;
              }
            }
          }
        }
      }
      else
      {
        // =======================
        // 3) In-situ statistics
        // =======================
        if (snapshot > 0 && indexTimeSample != 0 &&
            indexTimeSample % snapshot == 0)
        {
          auto startStats = system_clock::now();

          int nbNodes = m_mesh->getNumberOfNodes();
          for (int n = 0; n < nbNodes; ++n)
          {
            float v = pnGlobal(n, i1);
            minVal = std::min(minVal, v);
            maxVal = std::max(maxVal, v);
          }
        }
        auto endStats = system_clock::now();
      }
      auto endSnapshot = system_clock::now();
      snapshotTime_us +=
          duration_cast<microseconds>(endSnapshot - startSnapshot).count();
    }

    // =================================
    //  SNAPSHOT QUANTIZATION (ADD-ON) 
    // ===============================
    if (snapshot > 0 && indexTimeSample != 0 &&
    indexTimeSample % snapshot == 0 && !ad_hoc && quantification)
      {
          int nbNodes = m_mesh->getNumberOfNodes();
      
          QuantizationMeta qmeta;
          qmeta.minVal = minVal;
          qmeta.maxVal = maxVal;
          qmeta.nBits = nBit;
          
          
          std::vector<float> original(nbNodes);
          for (int n = 0; n < nbNodes; ++n)
              original[n] = pnGlobal(n, i1);
      
          // Quantification
          size_t elemSize = (nBit + 7)/8;   
          std::vector<uint8_t> buffer(nbNodes * elemSize);
      
          for (int n = 0; n < nbNodes; ++n)
          {
              float value = pnGlobal(n, i1);
              uint32_t q = quantize(value, qmeta);
              //if(value >= 0.05)
              //std::cout << q  << " "<< value << std::endl;
              std::memcpy(&buffer[n * elemSize], &q, elemSize);
          }
        
          // Reconstruction pour RMSE
          std::vector<float> reconstructed(nbNodes);
          for (int n = 0; n < nbNodes; ++n)
          {
              uint32_t q = 0;
              std::memcpy(&q, &buffer[n * elemSize], elemSize);
              reconstructed[n] = dequantize(q, qmeta);
              //std::cout << reconstructed[n] << std::endl;
          }
        
          // Calcul de la RMSE
          float rmse = computeRMSE(original, reconstructed);
          std::cout << "RMSE after quantization: " << rmse << std::endl;
        
          // Écrire dans le fichier binaire
          startComputeTime = system_clock::now();
          std::ofstream outQuant("snapshot_quant.bin", std::ios::binary | std::ios::app);
          outQuant.write(reinterpret_cast<char*>(&qmeta), sizeof(qmeta));
          outQuant.write(reinterpret_cast<char*>(buffer.data()), buffer.size());
          //std::cout << buffer.size() << std::endl;
          totalOutputTime += system_clock::now() - startOutputTime;

          originalSize += nbNodes * sizeof(float);
          compressedSize += buffer.size() * sizeof(uint8_t) + sizeof(qmeta);

      }


    // ==========================================
    //  Run-Length Encoding (RLE) for snapshots
    // ==========================================
    if (snapshot > 0 && indexTimeSample != 0 &&
        indexTimeSample % snapshot == 0 && !ad_hoc && RLE)
    {
      int nbNodes = m_mesh->getNumberOfNodes();
      QuantizationMeta qmeta{minVal, maxVal, nBit};

      //std::vector<uint16_t> quantizedValues(nbNodes);
      size_t elemSize = (nBit + 7)/8;   
      std::vector<uint8_t> buffer(nbNodes * elemSize);
      for (int n = 0; n < nbNodes; ++n)
      {
        float value = pnGlobal(n, i1);
        uint32_t q = quantize(value, qmeta);
        std::memcpy(&buffer[n * elemSize], &q, elemSize);
      }

      // Quantification déjà faite : quantizedValues
      auto rleData = RLE_compress(buffer,elemSize);
      startComputeTime = system_clock::now();
      writeRLE("snapshot_quant_rle.bin", rleData, qmeta);
      totalOutputTime += system_clock::now() - startOutputTime;


      // --- Mesure du taux de compression ---
      originalSize += nbNodes * sizeof(float);
      compressedSize += rleData.size() * sizeof(RLEPair) + sizeof(qmeta);
      
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

        if (!ad_hoc)
        {
          minPreRecv[m] = std::min(varnp1, minPreRecv[m]);
          maxPreRecv[m] = std::max(varnp1, maxPreRecv[m]);
          if (indexTimeSample == num_sample_ - 1)
          {
            meanRecv[m] += varnp1;
            meanRecv[m] /= static_cast<float>(num_sample_);
          }
          else
          {
            meanRecv[m] += varnp1;
          }
        }

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
  }
  // ============================
  // 4) Writting sismo in files
  // ============================
  if (recv_on) {


    if(ad_hoc){ //ad_hoc
      
      for (int m = 0; m < num_receivers; ++m) {

        const float xn = m_mesh->nodeCoord(rhsElementRcv[m], 0);
        const float yn = m_mesh->nodeCoord(rhsElementRcv[m], 1);
        const float zn = m_mesh->nodeCoord(rhsElementRcv[m], 2);

        const std::string filename = "recev_results_" + std::to_string(m) + ".csv";
        startComputeTime = system_clock::now();
        std::ofstream recv(filename, std::ios::app);

        if (!recv.is_open()) {
          std::cerr << "Error : opening this file is not possible" << filename << std::endl;
          continue;
        }

        
        if (recv.tellp() == 0) {
          recv << "indexTimeSample,xn,yn,zn,varnp1\n";
        }
        
        for (int i = 0; i < num_sample_; ++i) {
          recv << i << ","
               << xn << ","
               << yn << ","
               << zn << ","
               << pnAtReceiver(m, i) << "\n";
        }
        totalOutputTime += system_clock::now() - startOutputTime;

      }
    }

    if(RLE || quantification){
      

      size_t elemSize = (nBit + 7)/8;   
      
      for (int m = 0; m < num_receivers; ++m) {
        QuantizationMeta qmeta;
        qmeta.minVal = minPreRecv[m];
        qmeta.maxVal = maxPreRecv[m];
        qmeta.nBits = nBit;
        std::vector<uint8_t> buffer(num_sample_ * elemSize);
        std::vector<float> original(num_sample_);
        for (int i = 0; i < num_sample_; ++i) {
          float value = pnAtReceiver(m, i);
          uint32_t q = quantize(value, qmeta);
          std::memcpy(&buffer[i * elemSize], &q, elemSize);  
          original[i] = value;     
        }
        if(quantification){// quantification
          
               // Reconstruction pour RMSE
          std::vector<float> reconstructed(num_sample_);
          for (int n = 0; n < num_sample_; ++n)
          {
              uint32_t q = 0;
              std::memcpy(&q, &buffer[n * elemSize], elemSize);
              reconstructed[n] = dequantize(q, qmeta);
              //std::cout << reconstructed[n] << std::endl;
          }
        
          // Calcul de la RMSE
          float rmse = computeRMSE(original, reconstructed);
          startComputeTime = system_clock::now();
          std::cout << "RMSE after quantization: " << rmse << std::endl;
        
        // Écrire dans le fichier binaire
          std::ofstream outQuant("sismo_quant_" + std::to_string(m) + ".bin", std::ios::binary | std::ios::app);
          outQuant.write(reinterpret_cast<char*>(&qmeta), sizeof(qmeta));
          outQuant.write(reinterpret_cast<char*>(buffer.data()), buffer.size());
          totalOutputTime += system_clock::now() - startOutputTime;
          originalSize += num_sample_ * sizeof(float);
          compressedSize += buffer.size() * sizeof(uint8_t) + sizeof(qmeta);
        }
        if(RLE){
          startComputeTime = system_clock::now();
          auto rleData = RLE_compress(buffer,elemSize);
          writeRLE("sismo_quant_rle_"+ std::to_string(m) + ".bin", rleData, qmeta);
          totalOutputTime += system_clock::now() - startOutputTime;
          originalSize += num_sample_ * sizeof(float);
          compressedSize += rleData.size() * sizeof(RLEPair) + sizeof(qmeta);
        }
      }
    }
    
  } 
  auto totalEnd = system_clock::now();
  double totalRun_us =
      duration_cast<microseconds>(totalEnd - totalStart).count();

  // Convert to seconds
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
  if(!ad_hoc){
    if(snapshot){
      std::cout << "----   - Min by in-situ     : " << minVal << " " << sizeof(float) << " s\n";
      std::cout << "----   - Max by in-situ        : " << maxVal << " s\n";
      std::cout << "------------------------------------------------ \n";
    }
    if(RLE || quantification){
      double compressionRatio = static_cast<double>(originalSize) / compressedSize;
      std::cout << "---- Compression ratio: " << compressionRatio << "\n";
    }
    if(recv_on){
    for (int i = 0; i < num_receivers; i++)
      {
        std::cout << "------------------------------------------------ \n";
        std::cout << "----   - Min pressure of receiver " << i
                  << " by in-situ     : " << minPreRecv[i] << " s\n";
        std::cout << "----   - Max pressure of receiver " << i
                  << " by in-situ     : " << maxPreRecv[i] << " s\n";
        std::cout << "----   - Mean pressure of receiver " << i
                  << " by in-situ     : " << meanRecv[i] << " s\n";
        std::cout << "------------------------------------------------ \n";
      }
    }
  }
  std::cout << "----   - Sismo Time        : " << sismo_s << " s\n";
  std::cout << "------------------------------------------------ \n";
  std::cout << "---- Total run() Time      : " << total_s << " s\n";
  std::cout << "------------------------------------------------ \n";

  // -----------------------------------------
  // SAVE RUN TIME INFORMATION TO CSV
  // -----------------------------------------
  std::ofstream tfile("time_debug.csv", std::ios::app);
  if (!tfile)
  {
    std::cerr << "Error: cannot open time_debug.csv\n";
  }
  else
  {
    // If file is empty, write CSV header
    if (tfile.tellp() == 0)
    {
      tfile << "ex,ey,ez,snapshot,kernel_s,output_s,snapshot_s,sismo_s,total_"
               "s\n";
    }

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
