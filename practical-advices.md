## Performance

[] Pakiet BenchmarkTools posiada makro ``@btime`` , które wielokrotnie uruchamia kod i wyswietla rowniez alokacje pamieci,
``@benchmark`` wyswietla bardziej szczegolowe informacje

[] Julia najpierw iteruje po wartosciach w kolumnach a dopiero potem po wierszach (tak umieszcza je w pamieci).

[] ``@code_warntype baz()`` - wyswietla typy wystepujace w przebiegu funkcji baz()

[]
```Julia
julia> @macroexpand x = @. A + b*c
:(x = (+).(A, (*).(b, c)))
```

[] deklarowanie typow w funkcji nie jest koniecznie potrzebne kiedy dobrze jest napisany przebieg funckji, tak zeby komplilator wiedzial jakich typow uzyc. Jezeli ``@code_warntype`` zwroci unie typow przy wywolaniu oznacza to ze jest cos nie tak.

[] Profiler - wskazuje dokladny przebieg funkcji z informacja jakie linijki byly najdluzej wykonywane

[] Globalne wartosci deklarowac jako const (np. macierze ktore beda przekazywane do funkcji)

[] Funkcje inlinowe ``@inline``

[] Generowana funkcja na poziomie kompilacji kodu (74)

[] Subnormal numbers - numery bardzo bliskie zeru zamieniane na zero

[] Kiedy macierz jest np. transponowana, to typ takiej macierzy moze byc ``Adjoint`` czyli odpowiednik macierzy wejsciowej, to taka macierz najlepiej iterowac po wierszach

[] Do inicjalizacji macierzy wykorzystywac funkcje ``fill()``

[] ``@inbounds`` - wynikowy kod nie sprawdza czy odwolanie do tablicy nie zostalo przekroczone (a przeciez w wielu miejsach kodu wiemy ze sie nie przekroczy)

[] Prealokacja macierzy oraz przekazywanie macierzy przez referencje do funkcji

[] Wykorzysywac wersje funkcji z ``!`` ktora dziala na oryginalnej wersji tablicy, bez tworzenia kopii

[] Do operacji grupowych na tablicach wykorzysywac wersje funkcji z kropka ``.``, np. ``A .+ B``, ``sin.(A)``

[] ``@view`` pozwala podzielic macierz na czesci bez tworzenia kopii, np. ``@view(x[:,i])``

[] ``@simd`` - wykorzystanie mozliwosci wykonywania operacji o podobnej strutkruze (waruneki: to co sie dzieje miedzy iteracjami musi byc niezalezne i w srodku nie moze byc wywolyania funkcji) - boucd checking jest domyslnie wylaczone. Jest tez biblioteka **SIMD.jl**

[] Iterowac po obiektach a nie po indeksach ``for a in A``. A jezeli wartosci indeksu jest wymagany do oblczen to skorzystac z funkcji ``eachindex()``, np. ``for i in eachindex(A)``

[] Wielowatkowosc ``export JULIA_NUP_THREADS = 4``. ``@threads``

## Metoda potegowa:
- im wektor startowy jest bliższy poszukiwanemu wektorowi własnemu, tym mniej iteracji trzeba będzie wykonać.
- byc moze normalizacja nie jest az tak potrzebna
- dodac obliczanie pozostalych wartosci wlasnych

## Metoda QR
- wykorzystac Householdera
