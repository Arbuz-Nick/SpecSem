# Обзор библиотек и программ с фильтрацией облаков точек

### Point Cloud Library (PCL)
В первую очередь PCL это автономный крупномасштабный открытый проект для обработки 2D / 3D изображений и облаков точек. Проект начал зарождаться в 2011 году: на [GitHub-е](https://github.com/PointCloudLibrary "GitHub-е") проекта начали появляться первые коммиты. 
Хотя проект подходит для решения множества задач, в частности, нас будет интересовать раздел облаков точек, а именно их фильтрации. В PCL достаточно много различных фильтров для облаков точек. Например можно встретить такие фильтры как Гауссов фильтр ([StatisticalOutlierRemoval](https://pointclouds.org/documentation/classpcl_1_1_statistical_outlier_removal.html "StatisticalOutlierRemoval")), Двусторонний фильтр ([BilateralFilter](https://pointclouds.org/documentation/classpcl_1_1_bilateral_filter.html "BilateralFilter")), фильтр по количеству соседних точек в заданном радиусе ([RadiusOutlierRemoval](https://pointclouds.org/documentation/classpcl_1_1_radius_outlier_removal.html "RadiusOutlierRemoval")).
Разберем применение данной библиотеки на примере
```
int main(){
class A, B = 1
}
``` 
