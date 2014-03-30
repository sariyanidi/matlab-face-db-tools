% Sadece profil agiz kesmek icin duzenlenmis betik

clc
clear

imge_dizini = '/media/TOSHIBA EXT/data/session01/multiview';
nokta_dizini = '../noktalar/ilk_noktalar';
cikti_dizini = '../tmp/agiz';

noktalar = dir([nokta_dizini '/*.txt']);

pozlar = dir([nokta_dizini '/001_01_01_*_00c.txt']);
isiklar = {'00','01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19'}; %dir([nokta_dizini '/151_01_01_010_*.txt']);
secili_isiklar =  {'10','08','08','07','05','05','05','04','06','06','09','09','08'};

pzl = [1:5 7:12 14 15];%6 ve 13 numarali pozlari geciyoruz
% pzl = 1:15;

nokta_sayisi = [6 7 7 8 7 6 6 6 7 7 6 6 6];

iw = 90;
ih = 100;
%Data = zeros(iw*ih, 249*20*9);
%sayac = 1;
%=========POZ==========
for p = 1 : length(pzl)
    pz = pzl(p);
    poz = pozlar(pz).name(11:13);
    pozd = [poz(1:2) '_' poz(3)];
    
    disp(['Poz: ' poz])
    
    %=========KISI==========
    kisiler = dir([nokta_dizini '/*_01_01_' poz '_00c.txt']);
    for k = 1 : length(kisiler)
        kisi = kisiler(k).name(1:3);
        
        disp(['   Kisi: ' kisi])
    
        nokta = load([nokta_dizini '/' kisi '_01_01_' poz '_00c.txt']);
        if length(nokta) ~= nokta_sayisi(p)
            fprintf('%s kisisinde %s pozunda nokta sayilari hatali\n', kisi, poz);
            continue
        end
        %=========ISIK==========
        for is = 1 : length(secili_isiklar)
            isik = secili_isiklar{is};
%             isik = '00';
            imge = imread([imge_dizini '/' kisi '/01/' pozd '/' kisi '_01_01_' poz '_' isik '.png']);
            nokta = load([nokta_dizini '/' kisi '_01_01_' poz '_00c.txt']);
            
            disp(['      Isik: ' isik])
            
            switch p
                % agizlarin tamami sol tarafa baksin
                case {1,2,11,12,13}
                    imge = flipdim(imge,2);
                    nokta = nokta_fliple(nokta, size(imge,2));
                
                %sadece profil yuzlerden agiz kes
                case {3, 4, 10}
                    continue;
            end
            
            switch p
                case {2,5,9}
                    agizIdx = 5;
                otherwise
                    agizIdx = 4;
            end
            
            agiz_noktasi = nokta(agizIdx,:);
            gen = 34;
            dikd = [ agiz_noktasi(1)-gen+gen/5 agiz_noktasi(2)-gen 2*gen 2*gen ];
            agiz = imcrop(imge, dikd);
            agiz = imresize(agiz, [20 20], 'lanczos3');
            
            if 3 == length(size(agiz))
                agiz = rgb2gray(agiz);
            end
%             imshow(agiz,[]);
            imwrite(agiz, [cikti_dizini '/' kisi '_' poz '_' isik '.png']);
        end
    end
end
