# ****************************************************************** #
#            AVS Mesh Compression Exploration Model (MCEM)           #
# ****************************************************************** #

This is the description document of avs mesh compression exploration model.


# 1.build 
    ./build.sh

# 2.folder
    |---dependencies
        |---cmake
        |---directx-headers
        |---directx-math
        |---directx-mesh
        |---draco
        |---hpm-HPM-15.0
        |---nanoflann
        |---patches
        |---sal
        |---stb
        |---uvatlas
    |---source
        |---app
        |---lib
            |---common
            |---decode
            |---encode

# 3.dependencies
  # 3.1 directx-headers
  git clone https://github.com/microsoft/DirectX-Headers.git

  # 3.2 directx-math
  git clone https://github.com/microsoft/DirectXMath.git

  # 3.3 directx-mesh
  git clone https://github.com/microsoft/DirectXMesh.git

  # 3.4 draco
  git clone https://github.com/google/draco.git draco

  # 3.5 hpm
  Download from avs ftp remote library.

  # 3.6 stb
  git clone https://github.com/nothings/stb.git stb

  # 3.7 uvaltas
  git clone https://github.com/microsoft/UVAtlas.git

# 4.encode/decode config

  Those marked * are required.

  # encode
  
  Main options:
  -h | -?                           show help.
  -iMesh <input mesh>              *input mesh file name.
  -iMap <input map>                 input texture map file name.
  -oMesh <output mesh>              output mesh bitstream.
  -oMap <output map>                output texture map bitstream.
  -outputPath <output bitstream>    output total bitstream.
  -recMesh <reconstructed mesh>     output reconstructed mesh file.
  -encodeMap <flag>                 encode texture map(1) or not(0). default = 0.
  -inputPosBitdepth                 input position bitdepth.
  -removeDuplicate                  disable in lossless mode.
  -targetTriangleRatio              target triangle ratio in decimation.
  -subdivIterCount                  iteration count in subdivision fitting.

---------------------draco parameters---------------------------
  -point_cloud                      forces the input to be encoded as a point cloud.
  -draco_qp <value>                 quantization bits for the position attribute, default=11.
  -draco_qt <value>                 quantization bits for the texture coordinate attribute, default=10.
  -draco_qn <value>                 quantization bits for the normal vector attribute, default=8.
  -draco_qg <value>                 quantization bits for any generic attribute, default=8.
  -draco_cl <value>                 compression level [0-10], most=10, least=0, default=7.
  --skip ATTRIBUTE_NAME             skip a given attribute (NORMAL, TEX_COORD, GENERIC)
  --metadata                        use metadata to encode extra information in mesh files.
  -preserve_polygons                encode polygon info as an attribute.
Use negative quantization values to skip the specified attribute
---------------------hpm parameters---------------------------
  -hpm_cfg <hpm config>             hpm encode config path.
  -hpm_q <qp>                       hpm qp.
  -hpm_w <width>                    texture map width.
  -hpm_h <height>                   texture map height.
  -hpm_yuv <yuv file>               texture map yuv path. 
     

  # decode

  Main options:
  -h | -?                   show help.
  -i <input>               *input bitstream path.
  -oMesh <output>           output mesh file path.
  -oYuv <output>            output texture map yuv file path.
  -oMap <output>            output texture map png file path.

  # example

(lossy-geometry-lossy-texture)
-iMesh H:\database\category1_static\bread\bread_qp16_qt12.obj -encodeMap 1 -iMap H:\database\category1_static\bread\bread.jpg -hpm_cfg H:\database\mcem\dependencies\hpm-HPM-15.0\cfg\encode_AI.cfg -inputPosBitdepth 16 -hpm_q 45 -hpm_w 2048 -hpm_h 2048 -targetTriangleRatio 0.01 -subdivIterCount 1 -removeDuplicate 1 -oMesh H:\database\mcem\output\bread_c1_mesh.bin -oMap H:\database\mcem\output\bread_c1_map.bin -recMesh H:\database\mcem\output\bread_enc_c1.obj -outputPath H:\database\mcem\output\bread_c1.bin >H:\database\mcem\output\bread_enc_c1.log

-i H:\database\mcem\output\bread_c1.bin -oMesh H:\database\mcem\output\bread_dec_c1.obj -oYuv H:\database\mcem\output\bread_dec_c1.yuv -oMap H:\database\mcem\output\bread_dec_c1.png >H:\database\mcem\output\bread_dec_c1.log

(lossless-geometry-lossy-texture)
-iMesh H:\database\category1_static\bread\bread_qp16_qt12.obj -encodeMap 1 -iMap H:\database\category1_static\bread\bread.jpg -hpm_cfg H:\database\mcem\dependencies\hpm-HPM-15.0\cfg\encode_AI.cfg -inputPosBitdepth 16 -hpm_q 45 -hpm_w 2048 -hpm_h 2048 -targetTriangleRatio 1.0 -subdivIterCount 0 -removeDuplicate 0 -oMesh H:\database\mcem\output\bread_c2_mesh.bin -oMap H:\database\mcem\output\bread_c2_map.bin -recMesh H:\database\mcem\output\bread_enc_c2.obj -outputPath H:\database\mcem\output\bread_c2.bin >H:\database\mcem\output\bread_enc_c2.log

-i H:\database\mcem\output\bread_c2.bin -oMesh H:\database\mcem\output\bread_dec_c2.obj -oYuv H:\database\mcem\output\bread_dec_c2.yuv -oMap H:\database\mcem\output\bread_dec_c2.png >H:\database\mcem\output\bread_dec_c2.log
