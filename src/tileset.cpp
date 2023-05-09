#ifdef _WIN32
#include <gdal/ogr_spatialref.h>
#include <gdal/ogrsf_frmts.h>
#endif

#include <cmath>
#include <vector>
#include <string>
#include <cstring>
#include <algorithm>

#include <Eigen/Dense>
#include "extern.h"

///////////////////////
static const double pi = std::acos(-1);
static const int WGS84_EPSG = 4326;


#ifdef _WIN32
/**
 * 将输入EPSG投影下的坐标值转换为WGS84经纬度下的坐标值。
 *
 * @param insrs 输入坐标系的EPSG代码。
 * @param val 输入的坐标值，以double数组的形式给出，数组长度为2。
 * @param path 存储GDAL数据文件路径的字符串。
 * @return 转换成功则返回true，否则返回false。
 */
extern "C" bool epsg_convert(int insrs, double* val, char* path) {
    // 初始化输出坐标系
    static OGRSpatialReference outRs;
    static bool initialized = false;
    if (!initialized) {
        outRs.importFromEPSG(WGS84_EPSG);
        initialized = true;
    }

    // 设置GDAL数据文件路径
    CPLSetConfigOption("GDAL_DATA", path);
    
    // 初始化输入坐标系
    OGRSpatialReference inRs;
    inRs.importFromEPSG(insrs);

    // 创建坐标转换对象
    std::unique_ptr<OGRCoordinateTransformation> poCT(
        OGRCreateCoordinateTransformation(&inRs, &outRs));

    // 进行坐标转换
    if (poCT) {
        if (poCT->Transform(1, val, val + 1)) {
            return true;
        }
    }
    return false;
}

/**
 * 将输入的wkt坐标转换为WGS84经纬度坐标系。
 *
 * @param wkt WKT格式坐标系字符串。
 * @param val 输入坐标的数组指针，函数执行成功后存储转换后的经纬度坐标。
 * @param path GDAL_DATA环境变量的路径，用于GDAL库的数据路径。
 * @return 转换是否成功，成功返回true，否则返回false。
 */
extern "C" bool wkt_convert(char* wkt, double* val, char* path) {
    // 初始化WGS84坐标系
    static OGRSpatialReference outRs;
    static bool initialized = false;
    if (!initialized) {
        outRs.importFromEPSG(WGS84_EPSG);
        initialized = true;
    }

    // 设置GDAL_DATA环境变量
    CPLSetConfigOption("GDAL_DATA", path);

    // 创建OGRSpatialReference对象，并根据输入的WKT格式字符串初始化
    OGRSpatialReference inRs;
    inRs.importFromWkt(&wkt);

    // 创建OGRCoordinateTransformation对象，用于实现坐标转换
    std::unique_ptr<OGRCoordinateTransformation> poCT(
        OGRCreateCoordinateTransformation(&inRs, &outRs));

    // 执行坐标转换
    if (poCT) {
        if (poCT->Transform(1, val, val + 1)) {
            return true;
        }
    }
    return false;
}

#else
extern "C" bool
epsg_convert(int insrs, double* val, char* path)
{
    return false;
}

extern "C" bool
wkt_convert(char* wkt, double* val, char* path)
{
    return false;
}

#endif

extern "C"
{
    // 将角度转换为弧度的比例因子
    static const double DEG2RAD = 0.01745329251994329576923690768489;

    // 将经纬度转换为米的比例因子
    static const double LATI_TO_METER = 6335439.327292462539571136329447903;

    // 将经纬度转换为米的比例因子
    static const double LONGTI_TO_METER = 6378137.0 * DEG2RAD * std::cos(30.0 * DEG2RAD);

    // 将角度转换为弧度
    double degree2rad(double val) {
        return val * DEG2RAD;
    }

    // 将纬度转换为米
    double lati_to_meter(double diff) {
        return diff * LATI_TO_METER;
    }

    // 将经度转换为米
    double longti_to_meter(double diff, double lati) {
        return diff * LONGTI_TO_METER;
    }

    // 将米转换为纬度
    double meter_to_lati(double m) {
        return m / LATI_TO_METER;
    }

    // 将米转换为经度
    double meter_to_longti(double m, double lati) {
        return m / (LONGTI_TO_METER * std::cos(lati));
    }
}

/**

将给定的经纬度和高度转换为ECEF坐标系下的位置，并返回变换矩阵。

@param radian_x 经度，以弧度为单位

@param radian_y 纬度，以弧度为单位

@param height_min 相对于参考椭球体的高度，以米为单位

@return 返回16个元素的矩阵，表示从原点到变换后的点的变换矩阵
*/
std::vector<double> transfrom_xyz(double radian_x, double radian_y, double height_min) {
    // 定义椭球体参数
    const double ellipsod_a = 40680631590769;
    const double ellipsod_b = 40680631590769;
    const double ellipsod_c = 40408299984661.4;

    // 计算正弦和余弦
    const double cos_x = std::cos(radian_x);
    const double sin_x = std::sin(radian_x);
    const double cos_y = std::cos(radian_y);
    const double sin_y = std::sin(radian_y);

    // 计算经纬高转换为ECEF坐标系
    const double x0 = ellipsod_a * cos_y * cos_x;
    const double y0 = ellipsod_b * cos_y * sin_x;
    const double z0 = ellipsod_c * sin_y;
    const double gamma = std::sqrt(x0 * x0 + y0 * y0 + z0 * z0);

    // 计算单位向量
    const Eigen::Vector3d p(x0, y0, z0);
    const Eigen::Vector3d east_mat = (-y0 * p + x0 * Eigen::Vector3d::UnitZ()).normalized();
    const Eigen::Vector3d north_mat = (east_mat.cross(p)).normalized();

    // 计算转换矩阵
    Eigen::Matrix4d matrix = Eigen::Matrix4d::Identity();
    matrix.block<3, 1>(0, 0) = east_mat;
    matrix.block<3, 1>(0, 1) = north_mat;
    matrix.block<3, 1>(0, 2) = p.normalized();
    matrix(0, 3) = height_min * cos_y * cos_x;
    matrix(1, 3) = height_min * cos_y * sin_x;
    matrix(2, 3) = height_min * sin_y;

    // 转换矩阵转换为vector返回
    return std::vector<double>(matrix.data(), matrix.data() + matrix.size());
}

// 以C++的方式导出C函数，供其他语言调用
extern "C" void transform_c(double center_x, double center_y, double height_min, double* ptr) {
    // 将经纬度转换成弧度制
    double radian_x = degree2rad(center_x);
    double radian_y = degree2rad(center_y);
    // 调用函数计算ECEF坐标系下的变换矩阵
    std::vector<double> v = transfrom_xyz(radian_x, radian_y, height_min);
    // 将矩阵数据拷贝到指定内存中
    std::memcpy(ptr, v.data(), v.size() * sizeof(double));
}

/**

@brief 将模型的包围盒写入为 3D Tiles 标准的 json 文件，同时写入模型数据的 b3dm 文件路径
@param[in] trans 模型的位置、姿态信息
@param[in] box 模型的包围盒信息
@param[in] geometricError 模型的几何误差
@param[in] b3dm_file 模型数据的 b3dm 文件路径
@param[in] json_file 3D Tiles json 文件路径
@return bool 写入操作是否成功
@details 将给定的位置、姿态信息转换成仿射变换矩阵，然后将包围盒信息、几何误差、b3dm 文件路径
和仿射变换矩阵写入 3D Tiles json 文件中，最后将 json 文件和 b3dm 文件写入到磁盘上。
*/
bool write_tileset_box(
    Transform* trans, Box& box, 
    double geometricError, 
    const char* b3dm_file, 
    const char* json_file) {

    // 定义变量，用于存储矩阵和 JSON 文本
    std::vector<double> matrix;
    std::string json_txt = "{\"asset\": {\
        \"version\": \"0.0\",\
        \"gltfUpAxis\": \"Y\"\
    },\
    \"geometricError\":";
    std::string trans_str = "";

    // 判断是否需要进行坐标转换，如果需要则生成矩阵字符串
    if (trans) {
        matrix = transfrom_xyz(degree2rad(trans->radian_x), degree2rad(trans->radian_y), trans->min_height);
        trans_str = "\"transform\": [";
        for (int i = 0; i < 15; ++i) {
            trans_str += std::to_string(matrix[i]);
            if (i != 14) {
                trans_str += ",";
            }
        }
        trans_str += "1],";
        json_txt += trans_str;
    }

    // 生成 JSON 字符串中的 boundingVolume 字段
    json_txt += "\"boundingVolume\": {\"box\": [";
    for (int i = 0; i < 11; ++i) {
        json_txt += std::to_string(box.matrix[i]);
        json_txt += ",";
    }
    json_txt += std::to_string(box.matrix[11]);
    json_txt += "]},";

    // 生成 JSON 字符串中的 content 字段
    char last_buf[512];
    sprintf(last_buf, "\"geometricError\": %f,\
    \"refine\": \"REPLACE\",\
    \"content\": {\
    \"uri\": \"%s\"}}}", geometricError, b3dm_file);
    json_txt += last_buf;

    // 将 JSON 字符串写入文件
    bool ret = write_file(json_file, json_txt.c_str(), json_txt.size());
    if (!ret) {
        LOG_E("write file %s fail", json_file);
    }
    return ret;
}


bool write_tileset_region(
    Transform* trans,
    Region& region,
    double geometricError,
    const char* b3dm_file,
    const char* json_file)
{
    // 生成JSON字符串
    std::ostringstream json_txt;
    json_txt << R"({"asset": {"version": "0.0","gltfUpAxis": "Y"},"geometricError":)" << geometricError;
    json_txt << R"(,"root": {)";

    // 添加变换矩阵
    if (trans) {
        auto matrix = transfrom_xyz(trans->radian_x, trans->radian_y, trans->min_height);
        json_txt << R"("transform": [)";
        for (const auto& m : matrix) {
            json_txt << m << ",";
        }
        json_txt << R"(1],)";
    }

    // 添加包围盒信息
    json_txt << R"("boundingVolume": {"region": [)";
    auto pRegion = reinterpret_cast<double*>(&region);
    for (int i = 0; i < 5; ++i) {
        json_txt << pRegion[i] << ",";
    }
    json_txt << pRegion[5];

    // 添加B3DM文件链接
    json_txt << R"(]},"geometricError": )" << geometricError << R"(,"refine": "REPLACE","content": {"uri": ")";
    json_txt << b3dm_file << R"("}}})";

    // 将JSON字符串写入文件
    bool ret = write_file(json_file, json_txt.str().c_str(), json_txt.str().size());
    if (!ret) {
        LOG_E("write file %s fail", json_file);
    }
    return ret;
}

//写入tileset.json
bool write_tileset(
    double radian_x, double radian_y, 
    double tile_w, double tile_h, 
    double height_min, double height_max,
    double geometricError,
    const char* filename, const char* full_path)
{
    // 定义椭球体的参数
    double ellipsod_a = 40680631590769;
    double ellipsod_b = 40680631590769;
    double ellipsod_c = 40408299984661.4;

    // 计算投影坐标
    const double pi = std::acos(-1);
    double xn = std::cos(radian_x) * std::cos(radian_y);
    double yn = std::sin(radian_x) * std::cos(radian_y);
    double zn = std::sin(radian_y);

    double x0 = ellipsod_a * xn;
    double y0 = ellipsod_b * yn;
    double z0 = ellipsod_c * zn;
    double gamma = std::sqrt(xn*x0 + yn*y0 + zn*z0);
    double px = x0 / gamma;
    double py = y0 / gamma;
    double pz = z0 / gamma;
    double dx = x0 * height_min;
    double dy = y0 * height_min;
    double dz = z0 * height_min;

    // 计算东向和北向矩阵
    std::vector<double> east_mat = {-y0,x0,0};
    std::vector<double> north_mat = {
        (y0*east_mat[2] - east_mat[1]*z0),
        (z0*east_mat[0] - east_mat[2]*x0),
        (x0*east_mat[1] - east_mat[0]*y0)
    };

    // 计算东向和北向矩阵的模长
    double east_normal = std::sqrt(
        east_mat[0]*east_mat[0] + 
        east_mat[1]*east_mat[1] + 
        east_mat[2]*east_mat[2]
        );
    double north_normal = std::sqrt(
        north_mat[0]*north_mat[0] + 
        north_mat[1]*north_mat[1] + 
        north_mat[2]*north_mat[2]
        );

    // 归一化矩阵并生成转换矩阵
    std::vector<double> matrix = {
        east_mat[0] / east_normal,
        east_mat[1] / east_normal,
        east_mat[2] / east_normal,
        0,
        north_mat[0] / north_normal,
        north_mat[1] / north_normal,
        north_mat[2] / north_normal,
        0,
        xn,
        yn,
        zn,
        0,
        px + dx,
        py + dy,
        pz + dz,
        1
    };

    // 定义覆盖的区域
    std::vector<double> region = {
        radian_x - meter_to_longti(tile_w / 2, radian_y),
        radian_y - meter_to_lati(tile_h / 2),
        radian_x + meter_to_longti(tile_w / 2, radian_y),
        radian_y + meter_to_lati(tile_h / 2),
        0,
        height_max
    };

    std::string json_txt = "{\"asset\": {\
        \"version\": \"0.0\",\
        \"gltfUpAxis\": \"Y\"\
    },\
    \"geometricError\":";
    json_txt += std::to_string(geometricError);
    json_txt += ",\"root\": {\
        \"transform\": [";
    for (int i = 0; i < 15 ; i++) {
        json_txt += std::to_string(matrix[i]);
        json_txt += ",";
    }
    json_txt += "1],\
    \"boundingVolume\": {\
    \"region\": [";
    for (int i = 0; i < 5 ; i++) {
        json_txt += std::to_string(region[i]);
        json_txt += ",";
    }
    json_txt += std::to_string(region[5]);

    char last_buf[512];
    sprintf(last_buf,"]},\"geometricError\": %f,\
        \"refine\": \"REPLACE\",\
        \"content\": {\
        \"uri\": \"%s\"}}}", geometricError, filename);

    json_txt += last_buf;

    bool ret = write_file(full_path, json_txt.data(), (unsigned long)json_txt.size());
    if (!ret) {
        LOG_E("write file %s fail", filename);
    }
    return ret;
}