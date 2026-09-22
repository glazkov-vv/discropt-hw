# Как запускать

Данные о тестах и порогах указываются в файле `tests.txt`

Для запуска тестов используется скрипт `build_and_run.sh`. Результат выводится на экран и сохраняется в файле `report.txt`

# Отчет о работе

## Попытка 1. Отжиг в два этапа.

Задача решается в два этапа, как рекомендовано в условии: сначала распределяем клиентов по машинам, затем на каждой машине решаем TSP.

Первый этап (85% времени) — отжиг по распределению между машинами. Длина маршрута машины оценивается жадным обходом. Целевая функция этапа — сумма таких оценок по всем машинам.


Действия отжига:

- перемещение клиента в машину одного из 20 ближайших к нему клиентов;
- обмен двух близких клиентов из разных машин.

Они затрагивают только одну машину, так что пересчитывать каждый раз надо немного.

Начальная температура — средний положительный прирост от случайного перемещения, охлаждение экспоненциальное.

**Второй этап (15% времени) — отжиг маршрута каждой машины.** Стартуем с жадного обхода, ходы — 2-op и Or-opt. Оставшееся время делится равномерно между непустыми машинами.

Решение проходит все простые пороги и два сложных:

```
On test case ./data/vrp_16_3_1 got length 285.50. Easy limit passed. Score 5
On test case ./data/vrp_26_8_1 got length 607.65. Hard limit passed. Score 7
On test case ./data/vrp_51_5_1 got length 545.68. Easy limit passed. Score 5
On test case ./data/vrp_101_10_1 got length 830.68. Easy limit passed. Score 5
On test case ./data/vrp_200_16_1 got length 1412.79. Easy limit passed. Score 5
On test case ./data/vrp_421_41_1 got length 1912.46. Hard limit passed. Score 7
```

Слегка перераспределив распределение времени по этапам, получилось вщзять третий сложный тест.

```

On test case ./data/vrp_16_3_1 got length 285.50. Easy limit passed. Score 5
On test case ./data/vrp_26_8_1 got length 607.65. Hard limit passed. Score 7
On test case ./data/vrp_51_5_1 got length 545.68. Easy limit passed. Score 5
On test case ./data/vrp_101_10_1 got length 829.76. Hard limit passed. Score 7
On test case ./data/vrp_200_16_1 got length 1464.25. Easy limit passed. Score 5
On test case ./data/vrp_421_41_1 got length 1902.20. Hard limit passed. Score 7
```
