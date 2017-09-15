cmake_minimum_required(VERSION 3.7)

set(QT_USE_QTOPENGL TRUE)
set(QT_DEFINITIONS QwtPlot QwtWidgets QwtSvg)

find_package(Qt4 REQUIRED)

include(${QT_USE_FILE})

include_directories(/usr/include/qt4/Qt)

add_definitions(${QT_DEFINITIONS})

set(QWT_SOURCES qwt_abstract_scale_draw.cpp
                qwt_clipper.cpp
                qwt_color_map.cpp
                qwt_column_symbol.cpp
                qwt_date.cpp
                qwt_date_scale_draw.cpp
                qwt_date_scale_engine.cpp
                qwt_dyngrid_layout.cpp
                qwt_event_pattern.cpp
                qwt_graphic.cpp
                qwt_interval.cpp
                qwt_interval_symbol.cpp
                qwt_math.cpp
                qwt_magnifier.cpp
                qwt_null_paintdevice.cpp
                qwt_painter.cpp
                qwt_painter_command.cpp
                qwt_panner.cpp
                qwt_picker.cpp
                qwt_picker_machine.cpp
                qwt_pixel_matrix.cpp
                qwt_point_3d.cpp
                qwt_point_polar.cpp
                qwt_round_scale_draw.cpp
                qwt_scale_div.cpp
                qwt_scale_draw.cpp
                qwt_scale_map.cpp
                qwt_spline.cpp
                qwt_scale_engine.cpp
                qwt_symbol.cpp
                qwt_system_clock.cpp
                qwt_text_engine.cpp
                qwt_text_label.cpp
                qwt_text.cpp
                qwt_transform.cpp
                qwt_widget_overlay.cpp
                qwt_curve_fitter.cpp
                qwt_abstract_legend.cpp
                qwt_legend.cpp
                qwt_legend_data.cpp
                qwt_legend_label.cpp
                qwt_plot.cpp
                qwt_plot_renderer.cpp
                qwt_plot_xml.cpp
                qwt_plot_axis.cpp
                qwt_plot_curve.cpp
                qwt_plot_dict.cpp
                qwt_plot_directpainter.cpp
                qwt_plot_grid.cpp
                qwt_plot_histogram.cpp
                qwt_plot_item.cpp
                qwt_plot_abstract_barchart.cpp
                qwt_plot_barchart.cpp
                qwt_plot_multi_barchart.cpp
                qwt_plot_intervalcurve.cpp
                qwt_plot_zoneitem.cpp
                qwt_plot_tradingcurve.cpp
                qwt_plot_spectrogram.cpp
                qwt_plot_spectrocurve.cpp
                qwt_plot_scaleitem.cpp
                qwt_plot_legenditem.cpp
                qwt_plot_seriesitem.cpp
                qwt_plot_shapeitem.cpp
                qwt_plot_marker.cpp
                qwt_plot_textlabel.cpp
                qwt_plot_layout.cpp
                qwt_plot_canvas.cpp
                qwt_plot_panner.cpp
                qwt_plot_rasteritem.cpp
                qwt_plot_picker.cpp
                qwt_plot_zoomer.cpp
                qwt_plot_magnifier.cpp
                qwt_plot_rescaler.cpp
                qwt_point_mapper.cpp
                qwt_raster_data.cpp
                qwt_matrix_raster_data.cpp
                qwt_sampling_thread.cpp
                qwt_series_data.cpp
                qwt_point_data.cpp
                qwt_scale_widget.cpp
                qwt_plot_svgitem.cpp)

set(QWT_HEADERS qwt.h
                qwt_abstract_scale_draw.h
                qwt_clipper.h
                qwt_color_map.h
                qwt_compat.h
                qwt_column_symbol.h
                qwt_date.h
                qwt_date_scale_draw.h
                qwt_date_scale_engine.h
                qwt_dyngrid_layout.h
                qwt_global.h
                qwt_graphic.h
                qwt_interval.h
                qwt_interval_symbol.h
                qwt_math.h
                qwt_magnifier.h
                qwt_null_paintdevice.h
                qwt_painter.h
                qwt_painter_command.h
                qwt_panner.h
                qwt_picker.h
                qwt_picker_machine.h
                qwt_pixel_matrix.h
                qwt_point_3d.h
                qwt_point_polar.h
                qwt_round_scale_draw.h
                qwt_scale_div.h
                qwt_scale_draw.h
                qwt_scale_engine.h
                qwt_scale_map.h
                qwt_spline.h
                qwt_symbol.h
                qwt_system_clock.h
                qwt_text_engine.h
                qwt_text_label.h
                qwt_text.h
                qwt_transform.h
                qwt_widget_overlay.h
                qwt_curve_fitter.h
                qwt_event_pattern.h
                qwt_abstract_legend.h
                qwt_legend.h
                qwt_legend_data.h
                qwt_legend_label.h
                qwt_plot.h
                qwt_plot_renderer.h
                qwt_plot_curve.h
                qwt_plot_dict.h
                qwt_plot_directpainter.h
                qwt_plot_grid.h
                qwt_plot_histogram.h
                qwt_plot_item.h
                qwt_plot_abstract_barchart.h
                qwt_plot_barchart.h
                qwt_plot_multi_barchart.h
                qwt_plot_intervalcurve.h
                qwt_plot_tradingcurve.h
                qwt_plot_layout.h
                qwt_plot_marker.h
                qwt_plot_zoneitem.h
                qwt_plot_textlabel.h
                qwt_plot_rasteritem.h
                qwt_plot_spectrogram.h
                qwt_plot_spectrocurve.h
                qwt_plot_scaleitem.h
                qwt_plot_legenditem.h
                qwt_plot_seriesitem.h
                qwt_plot_shapeitem.h
                qwt_plot_canvas.h
                qwt_plot_panner.h
                qwt_plot_picker.h
                qwt_plot_zoomer.h
                qwt_plot_magnifier.h
                qwt_plot_rescaler.h
                qwt_point_mapper.h
                qwt_raster_data.h
                qwt_matrix_raster_data.h
                qwt_sampling_thread.h
                qwt_samples.h
                qwt_series_data.h
                qwt_series_store.h
                qwt_point_data.h
                qwt_scale_widget.h
                qwt_abstract_slider.h
                qwt_abstract_scale.h
                qwt_arrow_button.h
                qwt_analog_clock.h
                qwt_compass.h
                qwt_compass_rose.h
                qwt_counter.h
                qwt_dial.h
                qwt_dial_needle.h
                qwt_knob.h
                qwt_slider.h
                qwt_thermo.h
                qwt_wheel.h
                qwt_plot_svgitem.h)

add_library(${PROJECT_NAME}_plot STATIC
        ${QWT_HEADERS}
        ${QWT_SOURCES})

target_link_libraries(${PROJECT_NAME}_plot)

set(QWT_PLOT_SOURCES
        qwt_polar_fitter.cpp
        qwt_polar_item.cpp
        qwt_polar_picker.cpp
        qwt_polar_panner.cpp
        qwt_polar_magnifier.cpp
        qwt_polar_grid.cpp
        qwt_polar_curve.cpp
        qwt_polar_spectrogram.cpp
        qwt_polar_marker.cpp
        qwt_polar_itemdict.cpp
        qwt_polar_canvas.cpp
        qwt_polar_layout.cpp
        qwt_polar_renderer.cpp
        qwt_polar_plot.cpp)

set(QWT_POLT_HEADERS
        qwt_polar_global.h
        qwt_polar.h
        qwt_polar_fitter.h
        qwt_polar_item.h
        qwt_polar_picker.h
        qwt_polar_panner.h
        qwt_polar_magnifier.h
        qwt_polar_grid.h
        qwt_polar_curve.h
        qwt_polar_spectrogram.h
        qwt_polar_marker.h
        qwt_polar_itemdict.h
        qwt_polar_canvas.h
        qwt_polar_layout.h
        qwt_polar_renderer.h
        qwt_polar_plot.h
        )

add_library(${PROJECT_NAME}_plot_polar STATIC
        ${QWT_POLT_HEADERS}
        ${QWT_PLOT_SOURCES})

target_link_libraries(${PROJECT_NAME}_plot_polar)
























