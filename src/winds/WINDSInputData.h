/****************************************************************************
 * Copyright (c) 2025 University of Utah
 * Copyright (c) 2025 University of Minnesota Duluth
 *
 * Copyright (c) 2025 Behnam Bozorgmehr
 * Copyright (c) 2025 Jeremy A. Gibbs
 * Copyright (c) 2025 Fabien Margairaz
 * Copyright (c) 2025 Eric R. Pardyjak
 * Copyright (c) 2025 Zachary Patterson
 * Copyright (c) 2025 Rob Stoll
 * Copyright (c) 2025 Lucas Ulmer
 * Copyright (c) 2025 Pete Willemsen
 *
 * This file is part of QES-Winds
 *
 * GPL-3.0 License
 *
 * QES-Winds is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 *
 * QES-Winds is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QES-Winds. If not, see <https://www.gnu.org/licenses/>.
 ***************************************************************************/

/** @file WINDSInputData.h */

#pragma once
#include <cmath>

/*
 * A collection of data read from an XML. This contains
 * all root level information extracted from the xml.
 */

#include "util/ParseInterface.h"
#include "util/QESFileSystemHandler.h"

#include "SimulationParameters.h"
#include "FileOptions.h"
#include "MetParams.h"
#include "BuildingsParams.h"
#include "VegetationParams.h"

#include "fire/Fires.hpp"

#include "TURBParams.h"
#include "HRRRInput.h"

/**
 * @class WINDSInputData
 * @brief Collection of data read from an XML.
 *
 * Contains all root level information extracted from the xml
 *
 * @sa ParseInterface
 */
class WINDSInputData : public ParseInterface
{
public:
  SimulationParameters *simParams = nullptr; /**< :document this: */
  FileOptions *fileOptions = nullptr; /**< :document this: */
  MetParams *metParams = nullptr; /**< :document this: */
  TURBParams *turbParams = nullptr; /**< :document this: */
  BuildingsParams *buildingsParams = nullptr; /**< :document this: */
  VegetationParams *vegetationParams = nullptr; /**< :document this: */

  HRRRInput *hrrrInput = nullptr; /**< :HRRR input class instance: */

  Fires *fires = nullptr; /**< :document this: */

  WINDSInputData()
  {
    fileOptions = 0;
    metParams = 0;
    turbParams = 0;
    buildingsParams = 0;
    vegetationParams = 0;
    hrrrInput = 0;
    fires = 0;
  }

  WINDSInputData(const std::string fileName)
  {

    fileOptions = 0;
    metParams = 0;
    turbParams = 0;
    buildingsParams = 0;
    vegetationParams = 0;
    hrrrInput = 0;
    fires = 0;

    QESfs::set_file_path(fileName);
    // read and parse the XML
    parseXML(fileName, "QESWindsParameters");
  }

  /**
   * :document this:
   */
  virtual void parseValues()
  {
    parseElement<SimulationParameters>(true, simParams, "simulationParameters");
    parseElement<FileOptions>(false, fileOptions, "fileOptions");
    parseElement<MetParams>(false, metParams, "metParams");
    parseElement<TURBParams>(false, turbParams, "turbParams");
    parseElement<BuildingsParams>(false, buildingsParams, "buildingsParams");
    parseElement<VegetationParams>(false, vegetationParams, "vegetationParams");
    parseElement<HRRRInput>(false, hrrrInput, "hrrrInput");
    parseElement<Fires>(false, fires, "fires");
    if (metParams && metParams->windShapefile) {
      if (!metParams->sensors.empty() || !metParams->sensorName.empty()
          || simParams->wrfInputData || hrrrInput)
        throw ParseException("windShapefile cannot be combined with sensor, sensorName, WRF or HRRR input");
      if (simParams->totalTimeIncrements > 1)
        throw ParseException("windShapefile supports one timestep only");
      if (metParams->z0_domain_flag != 0 || simParams->domainRotation != 0
          || simParams->halo_x != 0 || simParams->halo_y != 0)
        throw ParseException("windShapefile requires uniform roughness, zero domainRotation and zero halo");
      if (simParams->originFlag != 0 && simParams->originFlag != 1)
        throw ParseException("windShapefile requires originFlag 0 or 1");
      for (int c = 0; c < 3; ++c)
        if (simParams->domain[c] < 2 || !std::isfinite(simParams->grid[c]) || simParams->grid[c] <= 0)
          throw ParseException("windShapefile requires at least two QES cells per axis and positive finite cellSize");
      if (simParams->verticalStretching != 0)
        throw ParseException("windShapefile currently supports uniform vertical QES grids only");
      if (buildingsParams && !buildingsParams->shpFile.empty()
          && (simParams->UTMx != static_cast<float>(metParams->windShapefile->originX)
              || simParams->UTMy != static_cast<float>(metParams->windShapefile->originY)))
        throw ParseException("windShapefile: building shapefile requires simulationParameters UTMx/UTMy matching originX/originY");
    }
  }
};
