---
title: "Tunnel diode 2024"
author: "YDS"
date: "2024-06-10"
output: html_document
---

# Введение

Туннельные диоды были открыты учёным Эсаки и используют эффект квантового туннелирования электронов через потенциальный барьер.

Этот эффект связан с волновой природой электронов, благодаря которой они могут попадать в классически запрещённые области.

Существует несколько видов туннельных диодов. Наиболее типичные сделаны на основе PN-переходов (рис. 1), но есть и другие (например, резонансно-туннельные диоды, диоды на основе свехрешёток и т.д.)

![Рисунок 1. Схема туннельного диода Эсаки [1].](images/esaki.png){width=400px}

$$ $$

Туннельные диоды известны тем, что на их ВАХ наблюдается область отрицательного дифференциального сопротивления (отрицательной дифференциальной проводимости --- ОДП). Дифференциальная проводимость определяется как:

$$ \Omega_{d}=\frac{dI}{dV} $$
Обычная проводимость --- отношение тока к напряжению --- не может быть отрицательной, а дифференциальная может. При этом ток положительный, но начинает падать при росте напряжения. Это можно видеть на рис. 2.

![Рисунок 2. ВАХ туннельного диода [1].](images/iv-curve.png){width=300px}

$$ $$

Участок ОДП может быть использован для генераторов переменного тока (генераторов колебаний), потому что он позволяет компенсировать внутреннее сопротивление цепи и избежать затухания колебаний.

# Теория

Ток через туннельный диод можно рассчитать по формуле Эсаки [2].

$$I(V)=\frac{em_{e}kT}{2\pi^{2}\hbar^{3}}\int_{0}^{\infty}T_{c}\left(E_{\perp}\right)\left[\ln\left(1+\exp\frac{E_{Fn}+eV-E_{\perp}}{kT}\right)-\ln\left(1+\exp\frac{E_{Fp}-E_{\perp}}{kT}\right)\right]\mathrm{d}E_{\perp}$$
Где $T_c$ --- вероятность (коэффициент) прохождения электронов через барьер, $E_{\perp}$ --- кинетическая энергия электронов в направлении границы между P и N областями, которая отсчитывается от дна зоны проводимости $E_c$.

Логарифмы описывают распределение электронов по энергиям в P и N областях.

Для работы туннельного диода области N и P должны быть сильно легированы, то там должно быть очень много примесей. Тогда при $V=0$ уровни Ферми будут находится внутри зоны проводимости в N области и внутри валентной зоны в P области, как изображено на рис. 3а.

![Рисунок 3. Зонная структура туннельного диода при а) $V=0$ и b) $V>0$ [3].](images/diode-scheme-2.png){width=500px}

$$ $$

Вероятность прохождения барьера может быть рассчитана с помощью квантовой механики для барьеров любой формы. В данном случае барьер похож на треугольный (рис. 4), но можно его считать и прямоугольным. При повышении напряжения уровни Ферми будут смещаться друг относительно друга и барьер между P и N областями будет уменьшаться (рис. 3b).

![Рисунок 4. Энергетический барьер.](images/diode-scheme-3.png){width=250px}

Самый простой вид имеет коэффициент прохождения для очень узкого барьера в виде дельта-функции:

$$U(x)=\alpha \delta(x),\qquad \alpha = Ha \tag{2}$$
где $H$ --- высота барьера в единицах энергии (эВ), а $a$ --- ширина барьера (например, в нм).

Для неё вероятность прохождения имеет вид:

$$T_c(E)=\frac{E}{E+\frac{m_e \alpha^2}{2\hbar^2}} \tag{3}$$


# Литература

1. Leo Esaki. Long Journey into Tunneling. Science, 22 March 1974, Volume 183, Number 4130.
1. N. Moulin, Mohamed Amara, F. Mandorlo, M. Lemiti. Tunnel junction I ( V ) characteristics: Review and a new model for p-n homojunctions. Journal of Applied Physics, 2019, 126 (3), pp.033105. 10.1063/1.5104314. hal-03035269
1. Messaadi Lotfi and Dibi Zohir. A Spice Behavioral Model of Tunnel Diode: Simulation and Application. International Journal of Control and Automation Vol. 9, No. 4 (2016), pp. 39-50 http://dx.doi.org/10.14257/ijca.2016.9.4.05
1. D. Mtn, M. PATIL, J. CHEN. Solid-State ElectronicsVol. 32, No. 1I, pp. 1025-1031, 1989