#!/usr/bin/python
# -*- coding: UTF-8 -*-

import sys
import numpy as np
from osgeo import gdal, osr
import math
import os
import time
import json
gdal.UseExceptions()

class Processor:
    water_level = 200

    def __init__(self, argv):
        if len(argv) < 2:
            print('用法: python dem_processor.py source_file')
            print('或输入 dem 文件路径')
            in_path = input('输入 dem 文件路径：')
        else:
            in_path = argv[1]
        self.timestamp = int(round(time.time() * 1000))
        self.readTif(in_path)
        # 淹没行列号
        self.flooded = []
        # 淹没数据标记，矩阵表示法，与 ReadAsArray 保持一致
        self.mark_data = np.zeros((self.dem_height, self.dem_width), dtype = int)

    def readTif(self, file_name):
        self.dirname = os.path.dirname(file_name)
        basename = os.path.basename(file_name)
        stem, suffix = os.path.splitext(basename)
        self.stem = stem
        self.suffix = suffix
        try:
            dataset = gdal.Open(file_name)
        except RuntimeError:
            print(f'没有找到文件：{file_name}')
            input('按任意键退出程序')
            sys.exit(1)
        if dataset == None:
            print(f'无法打开文件： {file_name}')
            input('按任意键退出程序')
            sys.exit(1)
        # 宽度，栅格矩阵的列数
        self.dem_width = dataset.RasterXSize
        # 高度，栅格矩阵的行数
        self.dem_height = dataset.RasterYSize
        # 波段数，蓝、绿、红、近红外
        self.dem_bands = dataset.RasterCount
        # 获取数据类型
        band1 = dataset.GetRasterBand(1)
        self.data_type = band1.DataType
        # 仿射矩阵信息
        self.dem_geotrans = dataset.GetGeoTransform()
        # 投影信息
        self.dem_proj = dataset.GetProjection()
        # 左上经度
        lt_lon = self.dem_geotrans[0]
        # 左上纬度
        lt_lat = self.dem_geotrans[3]
        # 右下经纬度
        rd_lon, rd_lat = self.tifRowcol2Lonlat(self.dem_width, self.dem_height)
        self.lt_lon = lt_lon
        self.lt_lat = lt_lat
        self.rd_lon = rd_lon
        self.rd_lat = rd_lat
        # 高程数据
        self.dem_data = dataset.ReadAsArray(0, 0, self.dem_width, self.dem_height)

    def tifRowcol2Lonlat(self, row, col):
        lon = self.dem_geotrans[0] + row * self.dem_geotrans[1] + col * self.dem_geotrans[2]
        lat = self.dem_geotrans[3] + row * self.dem_geotrans[4] + col * self.dem_geotrans[5]
        return lon, lat

    def tifLonlat2Rowcol(self, lon, lat):
        d_val = self.dem_geotrans[1] * self.dem_geotrans[5] - self.dem_geotrans[2] * self.dem_geotrans[4]
        row = (self.dem_geotrans[1] * (lat - self.dem_geotrans[3]) - self.dem_geotrans[4] * (lon - self.dem_geotrans[0])) / d_val
        col = (self.dem_geotrans[5] * (lon - self.dem_geotrans[0]) - self.dem_geotrans[2] * (lat - self.dem_geotrans[3])) / d_val
        row = math.floor(row)
        col = math.floor(col)
        return row, col

    def setStartPoint(self, lon, lat):
        self.start_lon = lon
        self.start_lat = lat
        row, col = self.tifLonlat2Rowcol(lon, lat)
        self.start_row = row
        self.start_col = col

    def reset(self):
        # 重置淹没行列号
        self.flooded = []
        # 重置淹没数据标记
        self.mark_data = np.zeros((self.dem_height, self.dem_width), dtype = int)

    def floodFill(self):
        '''
        八邻域填充算法
        '''
        self.flooded.append([self.start_row, self.start_col])
        flood_buffer_list = []
        flood_buffer_list.append([self.start_row, self.start_col])
        while len(flood_buffer_list) > 0:
            point = flood_buffer_list[0]
            point_row = point[0]
            point_col = point[1]
            flood_buffer_list.pop(0)
            # 向当前栅格单元的 8 个临近方向搜索连通域
            for r in [point_row - 1, point_row, point_row + 1]:
                for c in [point_col - 1, point_col, point_col + 1]:
                    row_index = r - 1
                    col_index = c - 1
                # 判断是否超出栅格边界，舍弃 tif 文件边界，防止某些 tif 不规范导致的异常 bug
                if (c >= self.dem_width or c <= 0 or r >= self.dem_height or r <= 0):
                    pass
                elif self.mark_data[row_index][col_index] == 1:
                    pass
                else:
                    current_elevation = self.dem_data[row_index][col_index]
                    # 如果低于水位
                    if current_elevation < self.water_level:
                        self.mark_data[row_index][col_index] = 1
                        self.flooded.append([r, c])
                        flood_buffer_list.append([r, c])

    def writeToFile(self):
        '''
        将淹没范围（边界点）以 GeoJSON 格式写入文件
        '''
        boundary_dict = {}
        boundary_point = []
        output_path = os.path.join(self.dirname, '.dem', f'{self.stem}_{self.water_level}.json')
        # 首次点位处理，将内部点去掉，仅保留边界点
        # 同时，为了防止存在线点（点连成后形成的是线段而非面），对点位做了一层延伸（上、左上、左三个邻域）
        first_res = []
        for item in self.flooded:
            row = item[0]
            col = item[1]
            row_index = row - 1
            col_index = col - 1
            add_flag = False
            if (row == self.dem_height or col == self.dem_width):
                add_flag = True
            elif (
                self.mark_data[row_index - 1][col_index] == 0 or
                self.mark_data[row_index + 1][col_index] == 0 or
                self.mark_data[row_index][col_index - 1] == 0 or
                self.mark_data[row_index][col_index + 1] == 0
            ):
                add_flag = True
            if add_flag:
                if (row > 1 and col > 1):
                    first_res.append([row - 1, col - 1])
                if row > 1:
                    first_res.append([row - 1, col])
                if col > 1:
                    first_res.append([row, col -1])
        # 将延伸点放入 flooded 中，并标记 mark_data
        for item in first_res:
            row = item[0]
            col = item[1]
            row_index = row - 1
            col_index = col - 1
            self.flooded.append([row, col])
            self.mark_data[row_index][col_index] = 1
        # 基于延伸后 flooded 进行点位处理，将内部点去掉，仅保留边界点
        for item in self.flooded:
            row = item[0]
            col = item[1]
            row_index = row - 1
            col_index = col - 1
            add_flag = False
            if (row == self.dem_height or col == self.dem_width):
                add_flag = True
            elif (
                self.mark_data[row_index - 1][col_index] == 0 or
                self.mark_data[row_index + 1][col_index] == 0 or
                self.mark_data[row_index][col_index - 1] == 0 or
                self.mark_data[row_index][col_index + 1] == 0
            ):
                add_flag = True
            if add_flag:
                boundary_dict[f'{row}_{col}'] = True
                boundary_point.append([row, col])
        if len(boundary_point) > 0:
            # 排序（如果这个淹没范围结果需要展示在页面上，那么 WebGIS 的加载需要点位的连续性，是顺时针或逆时针）
            # 选择第一个元素为起点，如果其周围有且只有两个点，则该点位生效，否则，选用下一个，直到找到符合的点
            # 这样做的目的是基于该点，顺时针或逆时针寻点，最终形成连续的点位
            start_index = 0
            start_point_effective = True
            while start_point_effective:
                start_point = boundary_point[start_index]
                start_row = start_point[0]
                start_col = start_point[1]
                side_point_list = []
                for r in [start_row - 1, start_row, start_row + 1]:
                    for c in [start_col - 1, start_col, start_col + 1]:
                        if (r == start_row and c == start_col):
                            pass
                        elif boundary_dict.get(f'{r}_{c}'):
                            side_point_list.append([r, c])
                if len(side_point_list) > 2:
                    start_index += 1
                else:
                    start_point_effective = False
            orderd_point = self.sortPoint(start_point, side_point_list, boundary_dict)
            orderd_boundary = []
            for item in orderd_point:
                row = item[0]
                col = item[1]
                lon, lat = self.tifRowcol2Lonlat(col, row)
                orderd_boundary.append([lon, lat])
            geojson = {
                "type": "FeatureCollection",
                "name": f'{self.water_level}_{len(orderd_boundary)}',
                "crs": {
                    "type": "name",
                    "properties": { "name": "urn:ogc:def:crs:EPSG::4490" }
                },
                "features": [
                    {
                        "type": "Feature",
                        "properties": { "NAME": " ", "REMARK": "" },
                        "geometry": {
                            "type": "Polygon",
                            "coordinates": [orderd_boundary]
                        }
                    }
                ]
            }
            with open(output_path, 'w') as output_file:
                output_file.write(json.dumps(geojson, indent=2))
            print(f'{self.water_level} 水位下的淹没数据量为 {len(self.flooded)} 边界数据量为 {len(orderd_boundary)}')

    def sortPoint(self, start_point, side_point_list, boundary_dict):
        second_point = side_point_list[1]
        orderd_point = []
        orderd_point.append(start_point)
        orderd_point.append(second_point)
        next_point = []
        # 开始寻点，直到形成闭环
        while next_point != side_point_list[0]:
            start_row = start_point[0]
            start_col = start_point[1]
            second_row = second_point[0]
            second_col = second_point[1]
            num = 1
            if (start_row == second_row + 1 and start_col == second_col - 1):
                num = 1
            elif (start_row == second_row + 1 and start_col == second_col):
                num = 2
            elif (start_row == second_row + 1 and start_col == second_col + 1):
                num = 3
            elif (start_row == second_row and start_col == second_col + 1):
                num = 4
            elif (start_row == second_row - 1 and start_col == second_col + 1):
                num = 5
            elif (start_row == second_row - 1 and start_col == second_col):
                num = 6
            elif (start_row == second_row - 1 and start_col == second_col - 1):
                num = 7
            elif (start_row == second_row and start_col == second_col - 1):
                num = 8
            next_point = self.getNextPoint(second_point, num, boundary_dict)
            orderd_point.append(next_point)
            start_point = second_point
            second_point = next_point
        return orderd_point

    def getNextPoint(self, center, start_num, boundary_dict):
        row = center[0]
        col = center[1]
        store = {
            'point1': [row + 1, col - 1],
            'point2': [row + 1, col],
            'point3': [row + 1, col + 1],
            'point4': [row, col + 1],
            'point5': [row - 1, col + 1],
            'point6': [row - 1, col],
            'point7': [row - 1, col - 1],
            'point8': [row, col - 1],
        }
        sum = 0
        point = []
        num = 0
        if start_num == 8:
            num = 1
        else:
            num = start_num + 1
        while sum != 8:
            item = store[f'point{num}']
            if num == 8:
                num = 1
            else:
                num += 1
            sum += 1
            if boundary_dict.get(f'{item[0]}_{item[1]}'):
                point = item
                sum = 8
        return point

    def run(self):
        self.reset()
        self.floodFill()
        self.writeToFile()

def main():
    process = Processor(sys.argv)
    # 假设当前水位为 306.35 ，输出 306.35 到 311.85 的淹没范围
    process.water_level = 306.35
    # 扩散点
    process.setStartPoint(118.1017278305, 36.1619407275)

    print(f'扩散起点: ({process.start_row}, {process.start_col})')
    start_elevation = process.dem_data[process.start_row - 1][process.start_col - 1]
    print(f'起点高程: {start_elevation}')

    while process.water_level <= 311.85:
        print(f'处理 {process.water_level} 水位的数据...')
        process.run()
        process.water_level = round(process.water_level + 0.01, 2)

    del process
    return 0

if __name__ == '__main__':
    sys.exit(main())
