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
 * ������EPSGͶӰ�µ�����ֵת��ΪWGS84��γ���µ�����ֵ��
 *
 * @param insrs ��������ϵ��EPSG���롣
 * @param val ���������ֵ����double�������ʽ���������鳤��Ϊ2��
 * @param path �洢GDAL�����ļ�·�����ַ�����
 * @return ת���ɹ��򷵻�true�����򷵻�false��
 */
extern "C" bool epsg_convert(int insrs, double* val, char* path) {
    // ��ʼ���������ϵ
    static OGRSpatialReference outRs;
    static bool initialized = false;
    if (!initialized) {
        outRs.importFromEPSG(WGS84_EPSG);
        initialized = true;
    }

    // ����GDAL�����ļ�·��
    CPLSetConfigOption("GDAL_DATA", path);
    
    // ��ʼ����������ϵ
    OGRSpatialReference inRs;
    inRs.importFromEPSG(insrs);

    // ��������ת������
    std::unique_ptr<OGRCoordinateTransformation> poCT(
        OGRCreateCoordinateTransformation(&inRs, &outRs));

    // ��������ת��
    if (poCT) {
        if (poCT->Transform(1, val, val + 1)) {
            return true;
        }
    }
    return false;
}

/**
 * �������wkt����ת��ΪWGS84��γ������ϵ��
 *
 * @param wkt WKT��ʽ����ϵ�ַ�����
 * @param val �������������ָ�룬����ִ�гɹ���洢ת����ľ�γ�����ꡣ
 * @param path GDAL_DATA����������·��������GDAL�������·����
 * @return ת���Ƿ�ɹ����ɹ�����true�����򷵻�false��
 */
extern "C" bool wkt_convert(char* wkt, double* val, char* path) {
    // ��ʼ��WGS84����ϵ
    static OGRSpatialReference outRs;
    static bool initialized = false;
    if (!initialized) {
        outRs.importFromEPSG(WGS84_EPSG);
        initialized = true;
    }

    // ����GDAL_DATA��������
    CPLSetConfigOption("GDAL_DATA", path);

    // ����OGRSpatialReference���󣬲����������WKT��ʽ�ַ�����ʼ��
    OGRSpatialReference inRs;
    inRs.importFromWkt(&wkt);

    // ����OGRCoordinateTransformation��������ʵ������ת��
    std::unique_ptr<OGRCoordinateTransformation> poCT(
        OGRCreateCoordinateTransformation(&inRs, &outRs));

    // ִ������ת��
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
    // ���Ƕ�ת��Ϊ���ȵı�������
    static const double DEG2RAD = 0.01745329251994329576923690768489;

    // ����γ��ת��Ϊ�׵ı�������
    static const double LATI_TO_METER = 6335439.327292462539571136329447903;

    // ����γ��ת��Ϊ�׵ı�������
    static const double LONGTI_TO_METER = 6378137.0 * DEG2RAD * std::cos(30.0 * DEG2RAD);

    // ���Ƕ�ת��Ϊ����
    double degree2rad(double val) {
        return val * DEG2RAD;
    }

    // ��γ��ת��Ϊ��
    double lati_to_meter(double diff) {
        return diff * LATI_TO_METER;
    }

    // ������ת��Ϊ��
    double longti_to_meter(double diff, double lati) {
        return diff * LONGTI_TO_METER;
    }

    // ����ת��Ϊγ��
    double meter_to_lati(double m) {
        return m / LATI_TO_METER;
    }

    // ����ת��Ϊ����
    double meter_to_longti(double m, double lati) {
        return m / (LONGTI_TO_METER * std::cos(lati));
    }
}

/**

�������ľ�γ�Ⱥ͸߶�ת��ΪECEF����ϵ�µ�λ�ã������ر任����

@param radian_x ���ȣ��Ի���Ϊ��λ

@param radian_y γ�ȣ��Ի���Ϊ��λ

@param height_min ����ڲο�������ĸ߶ȣ�����Ϊ��λ

@return ����16��Ԫ�صľ��󣬱�ʾ��ԭ�㵽�任��ĵ�ı任����
*/
std::vector<double> transfrom_xyz(double radian_x, double radian_y, double height_min) {
    // �������������
    const double ellipsod_a = 40680631590769;
    const double ellipsod_b = 40680631590769;
    const double ellipsod_c = 40408299984661.4;

    // �������Һ�����
    const double cos_x = std::cos(radian_x);
    const double sin_x = std::sin(radian_x);
    const double cos_y = std::cos(radian_y);
    const double sin_y = std::sin(radian_y);

    // ���㾭γ��ת��ΪECEF����ϵ
    const double x0 = ellipsod_a * cos_y * cos_x;
    const double y0 = ellipsod_b * cos_y * sin_x;
    const double z0 = ellipsod_c * sin_y;
    const double gamma = std::sqrt(x0 * x0 + y0 * y0 + z0 * z0);

    // ���㵥λ����
    const Eigen::Vector3d p(x0, y0, z0);
    const Eigen::Vector3d east_mat = (-y0 * p + x0 * Eigen::Vector3d::UnitZ()).normalized();
    const Eigen::Vector3d north_mat = (east_mat.cross(p)).normalized();

    // ����ת������
    Eigen::Matrix4d matrix = Eigen::Matrix4d::Identity();
    matrix.block<3, 1>(0, 0) = east_mat;
    matrix.block<3, 1>(0, 1) = north_mat;
    matrix.block<3, 1>(0, 2) = p.normalized();
    matrix(0, 3) = height_min * cos_y * cos_x;
    matrix(1, 3) = height_min * cos_y * sin_x;
    matrix(2, 3) = height_min * sin_y;

    // ת������ת��Ϊvector����
    return std::vector<double>(matrix.data(), matrix.data() + matrix.size());
}

// ��C++�ķ�ʽ����C���������������Ե���
extern "C" void transform_c(double center_x, double center_y, double height_min, double* ptr) {
    // ����γ��ת���ɻ�����
    double radian_x = degree2rad(center_x);
    double radian_y = degree2rad(center_y);
    // ���ú�������ECEF����ϵ�µı任����
    std::vector<double> v = transfrom_xyz(radian_x, radian_y, height_min);
    // ���������ݿ�����ָ���ڴ���
    std::memcpy(ptr, v.data(), v.size() * sizeof(double));
}

/**

@brief ��ģ�͵İ�Χ��д��Ϊ 3D Tiles ��׼�� json �ļ���ͬʱд��ģ�����ݵ� b3dm �ļ�·��
@param[in] trans ģ�͵�λ�á���̬��Ϣ
@param[in] box ģ�͵İ�Χ����Ϣ
@param[in] geometricError ģ�͵ļ������
@param[in] b3dm_file ģ�����ݵ� b3dm �ļ�·��
@param[in] json_file 3D Tiles json �ļ�·��
@return bool д������Ƿ�ɹ�
@details ��������λ�á���̬��Ϣת���ɷ���任����Ȼ�󽫰�Χ����Ϣ��������b3dm �ļ�·��
�ͷ���任����д�� 3D Tiles json �ļ��У���� json �ļ��� b3dm �ļ�д�뵽�����ϡ�
*/
bool write_tileset_box(
    Transform* trans, Box& box, 
    double geometricError, 
    const char* b3dm_file, 
    const char* json_file) {

    // ������������ڴ洢����� JSON �ı�
    std::vector<double> matrix;
    std::string json_txt = "{\"asset\": {\
        \"version\": \"0.0\",\
        \"gltfUpAxis\": \"Y\"\
    },\
    \"geometricError\":";
    std::string trans_str = "";

    // �ж��Ƿ���Ҫ��������ת���������Ҫ�����ɾ����ַ���
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

    // ���� JSON �ַ����е� boundingVolume �ֶ�
    json_txt += "\"boundingVolume\": {\"box\": [";
    for (int i = 0; i < 11; ++i) {
        json_txt += std::to_string(box.matrix[i]);
        json_txt += ",";
    }
    json_txt += std::to_string(box.matrix[11]);
    json_txt += "]},";

    // ���� JSON �ַ����е� content �ֶ�
    char last_buf[512];
    sprintf(last_buf, "\"geometricError\": %f,\
    \"refine\": \"REPLACE\",\
    \"content\": {\
    \"uri\": \"%s\"}}}", geometricError, b3dm_file);
    json_txt += last_buf;

    // �� JSON �ַ���д���ļ�
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
    // ����JSON�ַ���
    std::ostringstream json_txt;
    json_txt << R"({"asset": {"version": "0.0","gltfUpAxis": "Y"},"geometricError":)" << geometricError;
    json_txt << R"(,"root": {)";

    // ��ӱ任����
    if (trans) {
        auto matrix = transfrom_xyz(trans->radian_x, trans->radian_y, trans->min_height);
        json_txt << R"("transform": [)";
        for (const auto& m : matrix) {
            json_txt << m << ",";
        }
        json_txt << R"(1],)";
    }

    // ��Ӱ�Χ����Ϣ
    json_txt << R"("boundingVolume": {"region": [)";
    auto pRegion = reinterpret_cast<double*>(&region);
    for (int i = 0; i < 5; ++i) {
        json_txt << pRegion[i] << ",";
    }
    json_txt << pRegion[5];

    // ���B3DM�ļ�����
    json_txt << R"(]},"geometricError": )" << geometricError << R"(,"refine": "REPLACE","content": {"uri": ")";
    json_txt << b3dm_file << R"("}}})";

    // ��JSON�ַ���д���ļ�
    bool ret = write_file(json_file, json_txt.str().c_str(), json_txt.str().size());
    if (!ret) {
        LOG_E("write file %s fail", json_file);
    }
    return ret;
}

//д��tileset.json
bool write_tileset(
    double radian_x, double radian_y, 
    double tile_w, double tile_h, 
    double height_min, double height_max,
    double geometricError,
    const char* filename, const char* full_path)
{
    // ����������Ĳ���
    double ellipsod_a = 40680631590769;
    double ellipsod_b = 40680631590769;
    double ellipsod_c = 40408299984661.4;

    // ����ͶӰ����
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

    // ���㶫��ͱ������
    std::vector<double> east_mat = {-y0,x0,0};
    std::vector<double> north_mat = {
        (y0*east_mat[2] - east_mat[1]*z0),
        (z0*east_mat[0] - east_mat[2]*x0),
        (x0*east_mat[1] - east_mat[0]*y0)
    };

    // ���㶫��ͱ�������ģ��
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

    // ��һ����������ת������
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

    // ���帲�ǵ�����
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