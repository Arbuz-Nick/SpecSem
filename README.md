# Обзор библиотек и программ с фильтрацией облаков точек

Когда нам говорят "Облако точек", мы сразу можем представить что же это такое: это множество некоторых точек, у которых есть координаты, как-то сгруппированных в некоторое облако, которое что-то может из себя представлять: модель некоторого предмета или план местности.

![Пончик](https://upload.wikimedia.org/wikipedia/commons/thumb/4/4c/Point_cloud_torus.gif/220px-Point_cloud_torus.gif)

Чаще всего облака точек представляют собой только внешнюю поверхность модели, это связанно с методами их получения: облака точек создаются 3D-сканерами и фотограмметрическими методами обработки изображений. 3D-сканеры в автоматическом режиме замеряют большое количество точек на поверхности сканируемого объекта и зачастую генерируют на выходе облако точек в виде цифрового файла данных. 

Очевидно, что при работе любой программы, могут возникать различные ошибки, а когда в ней происходит огромное количество замеров, вероятность появления различный шумов и выбросов многократно увеличивается. 
Чтобы свести к минимуму количество таких ошибок, программисты написали огромное количество строк кода, создавая различные библиотеки и программные обеспечения для работы с облаками точек, их обработки и фильтрации. 
Ниже рассказывается о нескольких из них.

![Пример применения фильтра](https://russianblogs.com/images/294/f7cf2c1dd9e96414fc9ba5830a74d29e.png)


### Point Cloud Library (PCL)
![enter image description here](https://pointclouds.org/assets/images/logo.png)

В первую очередь PCL это автономный крупномасштабный открытый проект для обработки 2D / 3D изображений и облаков точек. Проект начал зарождаться в 2011 году: на [GitHub-е](https://github.com/PointCloudLibrary "GitHub-е") проекта начали появляться первые коммиты. 
Хотя проект подходит для решения множества задач, в частности, нас будет интересовать раздел облаков точек, а именно их фильтрации. В PCL достаточно много различных фильтров для облаков точек. Например можно встретить такие фильтры как Гауссов фильтр ([`StatisticalOutlierRemoval`](https://pointclouds.org/documentation/classpcl_1_1_statistical_outlier_removal.html "StatisticalOutlierRemoval")), Двусторонний фильтр ([`BilateralFilter`](https://pointclouds.org/documentation/classpcl_1_1_bilateral_filter.html "BilateralFilter")), фильтр по количеству соседних точек в заданном радиусе ([`RadiusOutlierRemoval`](https://pointclouds.org/documentation/classpcl_1_1_radius_outlier_removal.html "RadiusOutlierRemoval")). Все они объединяются в библиотеку **pcl_filters**, она содержит многие механизмы удаления выбросов и шума для трехмерных облаков точек.

Разберем применение данной библиотеки на примере двустороннего фильтра:

```c++
***
 pcl::PointCloud<pcl::PointXYZ>::Ptr xyz (new pcl::PointCloud<pcl::PointXYZ>); // Для представления облаков точек в библиотеке существует множество типов данных, например PointXYZ и PointCloud, о них будет ниже
 pcl::FastBilateralFilter<pcl::PointXYZ> fbf; // Создаем объект фильтра
 fbf.setInputCloud (xyz);
 fbf.setSigmaS (sigma_s); // Устанавливаем стандартное отклонение для пространственной окрестности двустороннего фильтра
 fbf.setSigmaR (sigma_r); // Устанавливаем стандартное отклонение по Гауссу, чтобы контролировать, насколько соседние пиксели падают из-за разницы в интенсивности
 pcl::PointCloud<pcl::PointXYZ> xyz_filtered; // Создаем объект, в который будем сохранять результат нашей фильтрации
 fbf.filter (xyz_filtered);
***
```
В первой строке мы упоминали тип данных [`PointXYZ`](https://pointclouds.org/documentation/structpcl_1_1_point_x_y_z.html), он и еще многие необходимые для работы типы (для представления 2D / 3D облаков, облаков с цветовой индикацией) лежат в библиотеке **pcl_common**.

Библиотека **pcl_common** содержит общие структуры данных и методы, используемые большинством библиотек PCL. Основные структуры данных включают класс PointCloud и множество типов точек, которые используются для представления точек, нормалей поверхности, значений цвета RGB, дескрипторов объектов и т. д. Она также содержит множество функций для вычисления расстояний / норм, средних и ковариаций, угловых преобразований, геометрических преобразований и многого другого. Помимо [`PointXYZ`](https://pointclouds.org/documentation/structpcl_1_1_point_x_y_z.html)(структуры точки с обычными евклидовыми координатами) в библиотеке есть такие типы точек, как [`PointXYZRGB`](https://pointclouds.org/documentation/structpcl_1_1_point_x_y_z_r_g_b.html) (евклидовы координаты + RGB цвет), [`PointXY`](https://pointclouds.org/documentation/structpcl_1_1_point_x_y.html)(точечная структура с двумя евклидовыми координатами), [`PointNormal`](https://pointclouds.org/documentation/structpcl_1_1_point_normal.html)(евклидовы координаты точки, нормальные координаты точки и оценка кривизны поверхности).

### Готовые приложения с функцией фильтрации облаков точек
Часто, при работе с облаками точек, нам нужно увидеть результат быстро и провзаимодействовать с ним, для этого уже написано несколько программ, частью функционала которых так же является и фильтрация облаков точек.

#### AutoCAD
![enter image description here](https://damassets.autodesk.net/content/dam/autodesk/www/product-imagery/lockup-610x66/autocad-lockup-610x66.png)

Одной из таких программ является AutoCAD, программа фирмы AutoDesk, являющейся крупнейшим поставщиком программного обеспечения для строительных и машиностроительных фирм. 
По своей сути и своему предназначению AutoCAD --- это программа для 2D или 3D проектирования и черчения. Первая версия была выпущена аж в 1982 году и обладала совсем небольшим числом элементарных объектов, таких как круги, дуги, линии и текст. Однако, на данный момент возможности AutoCAD значительно выросли, появилась и возможность работать с облаками точек, в том числе появилась возможность проводить их фильтрацию. 
В AutoCAD есть возможность фильтровать облака точек по:
- Классификации;
- Отметке;
- Значениям интенсивности LiDAR;
- Расположению на карте.


#### КРЕДО 3D СКАН
![enter image description here](https://credo-dialogue.ru/images/products/3D_skan/3dscan.png)

Данная отечественная программа предназначена для обработки облаков точек, полученных с использованием лазерного сканирования или фотограмметрическим методом и фотоизображений, полученных в процессе мобильного сканирования. Она предлагает довольно широкий спектр возможностей от простой работы с облаком точек до построения  цифровой модели рельефа, распознавания дорожных знаков, ЛЭП или уступов карьеров. 
Для фильтрации облаков в КРЕДО 3D СКАН существует несколько различных функций: 
- Прореживание облаков точек;
- Фильтр изолированных точек;
- Фильтр шумов ниже уровня рельефа.

**Прореживание облаков точек** нужно в тех случаях, когда на 3D модели присутствует некоторая плоская поверхность,  уменьшение плотности точек на которой не будет критичным, но позволит уменьшить затраты на обработку всей модели.
**Фильтр изолированных точек**, думаю, не нуждается в объяснении. Из модели просто убираются точки, которые в силу своей изолированности не играют особой роли или просто являются шумом или выбросом.
**Фильтр шумов ниже уровня рельефа** удобен конкретно для цели данной программы --- построения корректной цифровой модели рельефа
При его применении удаляются точки, которые лежат ниже рельефа (как это можно понять и из названия). Применение команды значительно ускоряет процесс выделения рельефа при будущей работе.

### Источники
- [Wiki](https://ru.wikipedia.org/wiki/%D0%9E%D0%B1%D0%BB%D0%B0%D0%BA%D0%BE_%D1%82%D0%BE%D1%87%D0%B5%D0%BA)
- [PCL](https://pointclouds.org/)
- [Russian Blogs](https://russianblogs.com/article/4262958562/)
- [AutoCAD](https://knowledge.autodesk.com/ru/support/revit/learn-explore/caas/CloudHelp/cloudhelp/2020/RUS/Revit-Model/files/GUID-BD499295-84DD-4BDE-B60D-73008AFBC791-htm.html)
- [КРЕДО 3D СКАН](https://credo-dialogue.ru/produkty/korobochnye-produkty/888-credo-3dscan-naznachenie.html)
