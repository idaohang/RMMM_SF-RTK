#ifndef MAINWINDOW_H
#define MAINWINDOW_H

//#include <QMainWindow>
#include <QGLWidget>
#include <QtOpenGL>
#include <mgl2/mgl.h>

class MainWindow : public QGLWidget
{
	Q_OBJECT
protected:
	mglGraph *gr;
	void resizeGL(int nWidth, int nHeight); // Метод вызываемый после каждого изменения размера окна
	void paintGL(); // Метод для вывода изображения на экран
void initializeGL(); // Метод для инициализирования opengl
public:
	MainWindow(QWidget *parent = 0);
	~MainWindow();
};

#endif // MAINWINDOW_H
