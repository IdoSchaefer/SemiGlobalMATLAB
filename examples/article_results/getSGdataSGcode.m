T = 1e3;
[allNt5, allmv5, aller5, max_ers5] = errorSGarticleSGcode(T, 5, 5, 4e3, 17);
pause(0.1)
[allNt7, allmv7, aller7, max_ers7] = errorSGarticleSGcode(T, 7, 7, 2.85e3, 14);
pause(0.1)
[allNt9, allmv9, aller9, max_ers9] = errorSGarticleSGcode(T, 9, 9, 2.85e3, 11);