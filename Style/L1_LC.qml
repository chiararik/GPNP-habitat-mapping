<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="3.28.12-Firenze" hasScaleBasedVisibilityFlag="0" minScale="1e+08" maxScale="0" styleCategories="AllStyleCategories">
  <flags>
    <Identifiable>1</Identifiable>
    <Removable>1</Removable>
    <Searchable>1</Searchable>
    <Private>0</Private>
  </flags>
  <temporal fetchMode="0" enabled="0" mode="0">
    <fixedRange>
      <start></start>
      <end></end>
    </fixedRange>
  </temporal>
  <elevation zscale="1" zoffset="0" enabled="0" symbology="Line" band="1">
    <data-defined-properties>
      <Option type="Map">
        <Option value="" type="QString" name="name"/>
        <Option name="properties"/>
        <Option value="collection" type="QString" name="type"/>
      </Option>
    </data-defined-properties>
    <profileLineSymbol>
      <symbol clip_to_extent="1" type="line" force_rhr="0" alpha="1" is_animated="0" frame_rate="10" name="">
        <data_defined_properties>
          <Option type="Map">
            <Option value="" type="QString" name="name"/>
            <Option name="properties"/>
            <Option value="collection" type="QString" name="type"/>
          </Option>
        </data_defined_properties>
        <layer class="SimpleLine" locked="0" pass="0" enabled="1">
          <Option type="Map">
            <Option value="0" type="QString" name="align_dash_pattern"/>
            <Option value="square" type="QString" name="capstyle"/>
            <Option value="5;2" type="QString" name="customdash"/>
            <Option value="3x:0,0,0,0,0,0" type="QString" name="customdash_map_unit_scale"/>
            <Option value="MM" type="QString" name="customdash_unit"/>
            <Option value="0" type="QString" name="dash_pattern_offset"/>
            <Option value="3x:0,0,0,0,0,0" type="QString" name="dash_pattern_offset_map_unit_scale"/>
            <Option value="MM" type="QString" name="dash_pattern_offset_unit"/>
            <Option value="0" type="QString" name="draw_inside_polygon"/>
            <Option value="bevel" type="QString" name="joinstyle"/>
            <Option value="183,72,75,255" type="QString" name="line_color"/>
            <Option value="solid" type="QString" name="line_style"/>
            <Option value="0.6" type="QString" name="line_width"/>
            <Option value="MM" type="QString" name="line_width_unit"/>
            <Option value="0" type="QString" name="offset"/>
            <Option value="3x:0,0,0,0,0,0" type="QString" name="offset_map_unit_scale"/>
            <Option value="MM" type="QString" name="offset_unit"/>
            <Option value="0" type="QString" name="ring_filter"/>
            <Option value="0" type="QString" name="trim_distance_end"/>
            <Option value="3x:0,0,0,0,0,0" type="QString" name="trim_distance_end_map_unit_scale"/>
            <Option value="MM" type="QString" name="trim_distance_end_unit"/>
            <Option value="0" type="QString" name="trim_distance_start"/>
            <Option value="3x:0,0,0,0,0,0" type="QString" name="trim_distance_start_map_unit_scale"/>
            <Option value="MM" type="QString" name="trim_distance_start_unit"/>
            <Option value="0" type="QString" name="tweak_dash_pattern_on_corners"/>
            <Option value="0" type="QString" name="use_custom_dash"/>
            <Option value="3x:0,0,0,0,0,0" type="QString" name="width_map_unit_scale"/>
          </Option>
          <data_defined_properties>
            <Option type="Map">
              <Option value="" type="QString" name="name"/>
              <Option name="properties"/>
              <Option value="collection" type="QString" name="type"/>
            </Option>
          </data_defined_properties>
        </layer>
      </symbol>
    </profileLineSymbol>
    <profileFillSymbol>
      <symbol clip_to_extent="1" type="fill" force_rhr="0" alpha="1" is_animated="0" frame_rate="10" name="">
        <data_defined_properties>
          <Option type="Map">
            <Option value="" type="QString" name="name"/>
            <Option name="properties"/>
            <Option value="collection" type="QString" name="type"/>
          </Option>
        </data_defined_properties>
        <layer class="SimpleFill" locked="0" pass="0" enabled="1">
          <Option type="Map">
            <Option value="3x:0,0,0,0,0,0" type="QString" name="border_width_map_unit_scale"/>
            <Option value="183,72,75,255" type="QString" name="color"/>
            <Option value="bevel" type="QString" name="joinstyle"/>
            <Option value="0,0" type="QString" name="offset"/>
            <Option value="3x:0,0,0,0,0,0" type="QString" name="offset_map_unit_scale"/>
            <Option value="MM" type="QString" name="offset_unit"/>
            <Option value="35,35,35,255" type="QString" name="outline_color"/>
            <Option value="no" type="QString" name="outline_style"/>
            <Option value="0.26" type="QString" name="outline_width"/>
            <Option value="MM" type="QString" name="outline_width_unit"/>
            <Option value="solid" type="QString" name="style"/>
          </Option>
          <data_defined_properties>
            <Option type="Map">
              <Option value="" type="QString" name="name"/>
              <Option name="properties"/>
              <Option value="collection" type="QString" name="type"/>
            </Option>
          </data_defined_properties>
        </layer>
      </symbol>
    </profileFillSymbol>
  </elevation>
  <customproperties>
    <Option type="Map">
      <Option value="false" type="bool" name="WMSBackgroundLayer"/>
      <Option value="false" type="bool" name="WMSPublishDataSourceUrl"/>
      <Option value="0" type="int" name="embeddedWidgets/count"/>
      <Option value="Value" type="QString" name="identify/format"/>
    </Option>
  </customproperties>
  <pipe-data-defined-properties>
    <Option type="Map">
      <Option value="" type="QString" name="name"/>
      <Option name="properties"/>
      <Option value="collection" type="QString" name="type"/>
    </Option>
  </pipe-data-defined-properties>
  <pipe>
    <provider>
      <resampling zoomedOutResamplingMethod="nearestNeighbour" maxOversampling="2" enabled="false" zoomedInResamplingMethod="nearestNeighbour"/>
    </provider>
    <rasterrenderer alphaBand="-1" nodataColor="" type="paletted" band="1" opacity="1">
      <rasterTransparency/>
      <minMaxOrigin>
        <limits>None</limits>
        <extent>WholeRaster</extent>
        <statAccuracy>Estimated</statAccuracy>
        <cumulativeCutLower>0.02</cumulativeCutLower>
        <cumulativeCutUpper>0.98</cumulativeCutUpper>
        <stdDevFactor>2</stdDevFactor>
      </minMaxOrigin>
      <colorPalette>
        <paletteEntry label="1 - Rocks, screes, debris" value="1" color="#bdbdbd" alpha="255"/>
        <paletteEntry label="2 - Snow and glaciers" value="2" color="#ccece6" alpha="255"/>
        <paletteEntry label="3 - Water" value="3" color="#084594" alpha="255"/>
        <paletteEntry label="4 - Broadleaved" value="4" color="#319812" alpha="255"/>
        <paletteEntry label="5 - Coniferous" value="5" color="#024618" alpha="255"/>
        <paletteEntry label="6 - Grassland" value="6" color="#97e800" alpha="255"/>
        <paletteEntry label="7 - Shrubs" value="7" color="#74c476" alpha="255"/>
        <paletteEntry label="8 - Wetland" value="8" color="#29c08e" alpha="255"/>
      </colorPalette>
      <colorramp type="randomcolors" name="[source]">
        <Option/>
      </colorramp>
    </rasterrenderer>
    <brightnesscontrast brightness="0" contrast="0" gamma="1"/>
    <huesaturation grayscaleMode="0" colorizeOn="0" colorizeGreen="128" colorizeBlue="128" colorizeStrength="100" saturation="0" colorizeRed="255" invertColors="0"/>
    <rasterresampler maxOversampling="2"/>
    <resamplingStage>resamplingFilter</resamplingStage>
  </pipe>
  <blendMode>0</blendMode>
</qgis>
