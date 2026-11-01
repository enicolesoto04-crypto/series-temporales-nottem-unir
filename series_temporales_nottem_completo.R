# LAB: SERIES TEMPORALES CON NOTTEM (1920-1939)
# Dataset incluido en R: no se necesita archivo CSV.
# Se trabaja con temperaturas mensuales en grados Fahrenheit.
# Ejecucion: abrir este archivo en RStudio y pulsar Source.

# Paquete necesario para la prueba de Dickey-Fuller aumentada.
if (!requireNamespace("tseries", quietly = TRUE)) {
  install.packages("tseries", repos = "https://cloud.r-project.org")
}
library(tseries)

# 1. CARGA Y EXPLORACION (20 %)
data("nottem", package = "datasets")
serie <- nottem
cat("\n=== 1. CARGA Y EXPLORACION ===\n")
cat("Clase:\n"); print(class(serie))
cat("Estructura:\n"); str(serie)
cat("Resumen:\n"); print(summary(serie))
cat("Inicio: "); print(start(serie))
cat("Fin: "); print(end(serie))
cat("Frecuencia: "); print(frequency(serie))
cat("Observaciones: ", length(serie), "\n", sep = "")
cat("Media: ", round(mean(serie), 2), " F\n", sep = "")
cat("Desviacion estandar: ", round(sd(as.numeric(serie)), 2), " F\n", sep = "")
cat("Minimo: ", min(serie), " F; maximo: ", max(serie), " F\n", sep = "")
plot(serie, col = "steelblue", lwd = 2,
     main = "Temperaturas mensuales de Nottingham (1920-1939)",
     xlab = "Ano", ylab = "Temperatura (grados Fahrenheit)")
grid()
# La forma repetida cada ano sugiere estacionalidad anual.

# 2. TENDENCIA, ESTACIONALIDAD Y RESIDUOS (30 %)
cat("\n=== 2. DESCOMPOSICION ===\n")
# Modelo aditivo: oscilacion estacional aproximadamente estable en grados F.
descomposicion <- decompose(serie, type = "additive")
plot(descomposicion)
cat("Efecto estacional promedio por mes (enero-diciembre):\n")
print(round(descomposicion$figure, 2))
# Componentes por separado con unidades y fechas originales.
plot(descomposicion$trend, main = "Tendencia suavizada", xlab = "Ano",
     ylab = "Temperatura (F)", col = "darkgreen", lwd = 2)
plot(descomposicion$seasonal, main = "Componente estacional anual", xlab = "Ano",
     ylab = "Desviacion estacional (F)", col = "darkorange", lwd = 2)
plot(descomposicion$random, main = "Componente irregular", xlab = "Ano",
     ylab = "Residuo (F)", col = "gray35")
cat("Interpretacion: los cambios mensuales siguen un ciclo anual de temperaturas\n")
cat("mas altas en meses calidos y mas bajas en meses frios; la tendencia\n")
cat("se examina en el componente trend sin asumir una subida sostenida.\n")

# 3. ESTACIONARIEDAD (30 %)
cat("\n=== 3. ESTACIONARIEDAD ===\n")
par(mfrow = c(1, 2))
acf(serie, main = "ACF - serie original")
pacf(serie, main = "PACF - serie original")
par(mfrow = c(1, 1))
adf_original <- adf.test(serie, alternative = "stationary")
cat("ADF de la serie original:\n"); print(adf_original)
cat("Hipotesis nula de ADF: raiz unitaria (no estacionariedad).\n")
cat("p < 0.05 aporta evidencia contra la hipotesis nula,\n")
cat("pero no descarta por si sola una estacionalidad marcada.\n")

# Si el ADF original no rechaza la raiz unitaria, aplicar diferencia simple.
if (adf_original$p.value >= 0.05) {
  serie_transformada <- diff(serie, differences = 1)
  transformacion <- "primera diferencia (lag = 1)"
} else {
  # Aun con ADF significativo, se elimina la periodicidad anual
  # para comparar el efecto de la estacionalidad sobre la ACF.
  serie_transformada <- diff(serie, lag = frequency(serie), differences = 1)
  transformacion <- "diferencia estacional (lag = 12)"
}
cat("Transformacion aplicada: ", transformacion, "\n", sep = "")
plot(serie_transformada, main = paste("Serie transformada:", transformacion),
     xlab = "Ano", ylab = "Cambio en temperatura (F)", col = "purple4")
par(mfrow = c(1, 2))
acf(serie_transformada, main = "ACF - serie transformada")
pacf(serie_transformada, main = "PACF - serie transformada")
par(mfrow = c(1, 1))
adf_transformada <- adf.test(serie_transformada, alternative = "stationary")
cat("ADF de la serie transformada:\n"); print(adf_transformada)
comparacion_adf <- data.frame(
  Serie = c("Original", "Transformada"),
  Estadistico_ADF = c(unname(adf_original$statistic), unname(adf_transformada$statistic)),
  p_valor = c(adf_original$p.value, adf_transformada$p.value))
cat("Tabla comparativa ADF (antes/despues):\n")
print(comparacion_adf, row.names = FALSE)
cat("Resultado despues de transformar: ",
    if (adf_transformada$p.value < 0.05)
      "ADF rechaza raiz unitaria al 5%.\n"
    else "ADF no rechaza raiz unitaria al 5%.\n", sep = "")
# Una prueba ADF no garantiza que la serie carezca de estacionalidad.

# 4. VALORES ATIPICOS (20 %)
cat("\n=== 4. DETECCION DE POSIBLES OUTLIERS ===\n")
boxplot(as.numeric(serie), horizontal = TRUE, col = "lightblue",
        main = "Boxplot de temperaturas originales", xlab = "Temperatura (F)")
# Importante: temperaturas altas en verano o bajas en invierno pueden ser normales.
# Por ello buscamos puntos atipicos en los residuos de la descomposicion.
residuos <- descomposicion$random
res_validos <- as.numeric(residuos[!is.na(residuos)])
q <- quantile(res_validos, probs = c(0.25, 0.75))
rango_iqr <- unname(q[2] - q[1])
limite_inferior <- unname(q[1] - 1.5 * rango_iqr)
limite_superior <- unname(q[2] + 1.5 * rango_iqr)
indices_atipicos <- which(!is.na(residuos) &
                          (residuos < limite_inferior | residuos > limite_superior))
cat("Limites IQR para residuos: ", round(limite_inferior, 2), " y ",
    round(limite_superior, 2), " F\n", sep = "")
cat("Numero de posibles atipicos en residuos: ", length(indices_atipicos), "\n", sep = "")
plot(serie, main = "Posibles valores atipicos ajustados por estacionalidad",
     xlab = "Ano", ylab = "Temperatura (F)", col = "steelblue", lwd = 1.5)
if (length(indices_atipicos) > 0) {
  points(time(serie)[indices_atipicos], as.numeric(serie)[indices_atipicos],
         col = "red", pch = 19, cex = 1.1)
  tabla_atipicos <- data.frame(
    Ano = floor(time(serie)[indices_atipicos]),
    Mes = cycle(serie)[indices_atipicos],
    Temperatura_F = as.numeric(serie)[indices_atipicos],
    Residuo_F = round(as.numeric(residuos)[indices_atipicos], 2))
  print(tabla_atipicos, row.names = FALSE)
} else {
  cat("No se detectaron observaciones fuera del umbral IQR residual.\n")
}
legend("topright", legend = c("Temperatura observada", "Posible atipico (IQR residual)"),
       col = c("steelblue", "red"), lty = c(1, NA), pch = c(NA, 19), bty = "n", cex = 0.8)
cat("Interpretacion: los puntos rojos son candidatos estadisticos;\n")
cat("no se puede afirmar que sean errores sin contexto meteorologico adicional.\n")

# 5. CONCLUSIONES
cat("\n=== 5. CONCLUSIONES ===\n")
cat("La serie contiene temperaturas mensuales de Nottingham de 1920 a 1939.\n")
cat("Se observa una estacionalidad anual; la descomposicion aditiva permite\n")
cat("separarla de la tendencia de largo plazo y de la variacion irregular.\n")
cat("La ACF/PACF y las pruebas ADF documentan el comportamiento original\n")
cat("y el obtenido tras la transformacion indicada.\n")
cat("Los atipicos se buscaron en los residuos, evitando confundir los\n")
cat("meses calidos o frios normales con anomalías.\n")
