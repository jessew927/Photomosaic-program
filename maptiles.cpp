/**
* @file maptiles.cpp
* Code for the maptiles function.
*/

#include <iostream>
#include <map>
#include "maptiles.h"
//#include "cs225/RGB_HSL.h"

using namespace std;

Point<3> convertToXYZ(LUVAPixel pixel) {
  return Point<3>( pixel.l, pixel.u, pixel.v );
}

MosaicCanvas* mapTiles(SourceImage const& theSource, vector<TileImage>& theTiles)
{

    MosaicCanvas* mosaic = new MosaicCanvas(theSource.getRows(),theSource.getColumns());
    vector<Point<3>> points;
    map<Point<3>, TileImage*> match;

    for (auto it = theTiles.begin(); it!=theTiles.end(); ++it) {
      points.push_back(convertToXYZ((*it).getAverageColor()));
      match.insert({ convertToXYZ((*it).getAverageColor()), &(*it)});
    }

    KDTree<3> images(points);

    for(int i = 0;i<theSource.getRows();i++){
      for(int j = 0;j<theSource.getColumns();j++){
        Point<3> replacement = images.findNearestNeighbor(convertToXYZ(theSource.getRegionColor(i,j)));
        size_t count = match.count(replacement);
        if(count > 0){
          mosaic->setTile(i,j,match[replacement]);
        }
      }
    }

    return mosaic;
  }
