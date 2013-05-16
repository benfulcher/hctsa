function s = macro(s) 
s = cut(s, 2, 2, 2); 
s = embed(s, 10, 7, 1, 'Rect'); 
s = corrdim2(s, 1000, 0.1, 20); 
