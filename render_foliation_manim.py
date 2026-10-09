#!/usr/bin/env python3
"""
Manim Scene: 17-Layer Modular Foliation Dynamics & Radion Wavefunction
Run with: manim -pql render_foliation_manim.py FoliationBulkScene
"""

from manim import *
import json
import numpy as np

class FoliationBulkScene(ThreeDScene):
    def construct(self):
        # 1. Configuración de cámara y entorno 3D
        self.set_camera_orientation(phi=70 * DEGREES, theta=-55 * DEGREES)
        
        title = Title("17-Layer Modular Bulk Foliation & Radion Localization", font_size=28)
        self.add_fixed_in_frame_mobjects(title)

        # Cargar datos numéricos
        with open("foliation_geometry_3d.json") as f:
            data = json.load(f)

        y_nodes = data["y_nodes"]
        n_layers = len(y_nodes) - 1

        # 2. Ejes del Bulk (y: dimensión extra, x,z: espacio 4D)
        axes = ThreeDAxes(
            x_range=[-0.5, 7.5, 1],
            y_range=[-2.5, 2.5, 1],
            z_range=[-2.5, 2.5, 1],
            x_length=8.0,
            y_length=4.5,
            z_length=4.5
        )
        self.play(Create(axes), run_time=1.5)

        # Etiquetas UV e IR
        uv_label = Text("UV Brane (Planck)", font_size=18, color=BLUE_B).move_to(axes.c2p(0, 0, 2.7))
        ir_label = Text("IR Brane (TeV)", font_size=18, color=RED_B).move_to(axes.c2p(7.0, 0, 2.7))
        self.add_fixed_in_frame_mobjects(uv_label, ir_label)

        # 3. Construcción de las 17 láminas (branas / interfaces)
        layers_group = VGroup()
        for idx in range(n_layers):
            y_coord = (y_nodes[idx] / 35.0) * 7.0
            t = idx / float(n_layers)
            col = interpolate_color(BLUE, RED, t)
            
            # Plano transversal
            plane = Square(
                side_length=3.2,
                fill_color=col,
                fill_opacity=0.18 + 0.15 * t,
                stroke_color=col,
                stroke_width=1.2
            )
            # Orientar el plano perpendicular al eje X (que representa 'y')
            plane.rotate(PI / 2, axis=UP)
            plane.move_to(axes.c2p(y_coord, 0, 0))
            layers_group.add(plane)

        self.play(LaggedStartMap(FadeIn, layers_group, lag_ratio=0.08), run_time=3.0)

        # 4. Función de onda del Radión (Modo escalar normalizado hacia la brana IR)
        radion_curve = ParametricFunction(
            lambda t: axes.c2p(
                t,
                0,
                1.8 * np.exp(0.4 * (t - 7.0)) * np.cos(3.0 * t)
            ),
            t_range=[0, 7.0],
            color=YELLOW,
            stroke_width=4
        )
        radion_tag = Text("Scalar Radion: 602.50 GeV", font_size=20, color=YELLOW).move_to(axes.c2p(4.5, 0, -2.2))
        self.add_fixed_in_frame_mobjects(radion_tag)

        self.play(Create(radion_curve), run_time=2.0)

        # 5. Rotación cinemática continua de cámara
        self.begin_ambient_camera_rotation(rate=0.20)
        self.wait(4.0)
        self.stop_ambient_camera_rotation()

        self.play(FadeOut(radion_curve), FadeOut(layers_group), run_time=1.5)
