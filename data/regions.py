import ee 
import geemap.core as geemap


'''
Regions to Define 
- Raleigh: -78.69649427302112, 35.7664191909458 (SOUTHWEST RALEIGH)
- Boone:  -81.68516978460896, 36.21834246106182 
- Asheville:  -82.55831963941482, 35.633546529345516 (NORTH AVL)
- Charleston: -79.9294956284804, 32.796133071541156
- Shennandoah National Park: -78.40685152551191, 38.499092552797826
'''

#buffer distance in meters
UNIV_BUFFER = 40000

Raleigh = ee.Geometry.Point([-78.69649427302112, 35.7664191909458]).buffer(UNIV_BUFFER)

Boone = ee.Geometry.Point([-81.68516978460896, 36.21834246106182]).buffer(UNIV_BUFFER)

Asheville = ee.Geometry.Point([-82.55831963941482, 35.633546529345516]).buffer(UNIV_BUFFER)

Charleston = ee.Geometry.Point([-79.9294956284804, 32.796133071541156]).buffer(UNIV_BUFFER)

Shennandoah = ee.Geometry.Point([-78.40685152551191, 38.499092552797826]).buffer(UNIV_BUFFER)

Map = geemap.Map()




