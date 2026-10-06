#!/usr/bin/env python3
"""Build the plain-language review page (easy_review_KR.html) with embedded figures."""
import base64
import os, os, html

ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
B = f'{ANALYSIS}/gz'
ASSETS = f'{B}/figures/review_assets'
FIG = f'{ASSETS}/easyfig'
OUT = f'{B}/results/review/easy_review_KR.html'


PUB = f'{B}/figures'


def refresh():
    """Re-encode every rendered figure, so the page can never show a figure that has since been redrawn."""
    from PIL import Image
    os.makedirs(FIG, exist_ok=True)
    for n in [f'Figure{i}' for i in range(1, 6)] + [f'FigureS{i}' for i in range(1, 8)]:
        src = f'{PUB}/{"supp" if "S" in n[6:] else "pub"}/{n}.png'
        key = n.replace('Figure', 'F').replace('FS', 'S')
        im = Image.open(src).convert('RGB')
        w = 1400
        im.resize((w, round(im.height * w / im.width)), Image.LANCZOS).save(
            f'{FIG}/{key}.jpg', 'JPEG', quality=82, optimize=True)
    print('figures re-encoded:', ', '.join(sorted(os.listdir(FIG))))


refresh()


def img(key):
    with open(f'{FIG}/{key}.jpg', 'rb') as fh:
        return 'data:image/jpeg;base64,' + base64.b64encode(fh.read()).decode()


# ---------------------------------------------------------------- content ----

MAIN_FIGS = [
    dict(
        key='F1', num='그림 1', role='있는가', rolecls='r1',
        q='하락이 정말 있는가, 그리고 언제 일어나는가?',
        what='네 종(생쥐·소·돼지·토끼)의 난자부터 상실배까지, 그리고 생쥐 두 번째 데이터에서 시계 값이 어떻게 움직이는지. '
             '마지막 칸은 “모계 RNA가 사라지는 것만으로 이 하락이 만들어지는가”를 시험한 결과.',
        why='논문의 출발점을 세우는 그림이다. 두 가지를 동시에 해야 한다. ① 하락이 한 연구의 우연이 아니라 독립된 '
            '생쥐 데이터에서 반복된다는 것, ② 그 하락 구간이 <b>세포분열이 없는</b> 2세포기 안이라는 것. '
            '분열이 없으므로 “세포가 바뀌어서 값이 바뀐 것”이라는 가장 흔한 반론이 여기서 미리 차단된다.',
        without='이 그림이 없으면 리뷰어는 “한 데이터셋의 아티팩트 아니냐”, “세포 조성이 바뀐 것 아니냐”를 가장 먼저 묻는다. '
                '특히 (d)가 없으면 논문 전체가 조성 착시 하나로 무너질 수 있다.',
        panels=[('a', '연구 설계 모식도 — 어느 단계에서 minor/major ZGA가 일어나고, 어느 데이터가 어느 단계를 덮는지. '
                      '초기 2세포기와 후기 2세포기 사이에 분열이 없다는 점을 그림으로 못박는다.'),
                ('b', '4종 시계 궤적 — 사전에 정한 첫 기준(“최대 하락이 문헌상 major ZGA 구간과 겹치는가”)의 결과. '
                      '<b>생쥐만 충족, 4종 중 1종.</b> 음성 결과를 숨기지 않고 본문 첫 그림에 그대로 둔 패널.'),
                ('c', '독립 재현(GSE45719) — 검은 점은 전체 유전자(V0), 파란 점은 동적 유전자를 뺀 값(V2). '
                      '유전자를 빼도 하락이 남는다는 것을 한 그림에서 보여준다.'),
                ('d', '조성 검정 — 모계 전사체를 관찰된 비율만큼 지우기만 해도 원래 하락의 0.70(생쥐)·0.47(소)이 재현된다. '
                      '즉 “절반쯤은 청소 효과”라고 논문이 먼저 인정하는 패널.')],
        nums=['생쥐 초기→후기 2세포기 −0.19 (최대 하락 구간일 확률 1.00)',
              'GSE45719 재현 −0.14 · GSE66582 −0.26',
              '사전 기준 충족 4종 중 1종', '모계 제거 모의 재현율 0.70 / 0.47'],
    ),
    dict(
        key='F2', num='그림 2', role='의존하는가', rolecls='r2',
        q='그 하락이 ZGA(접합자 유전체 활성화)에 의존하는가?',
        what='ZGA를 막은 배아와 막지 않은 배아를 나란히 놓고, 각 arm에서 초기→후기 2세포기 하락이 얼마나 되는지. '
             '라이브러리 하나가 점 하나다.',
        why='이 논문이 기존 연구와 갈라지는 지점. 지금까지의 배아 시계 연구는 모두 <b>정상 배아를 관찰</b>하기만 했다. '
            '여기서는 ZGA를 실제로 막았을 때 하락이 어떻게 되는지를 본다. '
            '세 가지 기호 D(각 arm의 하락), I(교란 − 대조), R(구제 − 교란)가 논문 전체에서 계속 쓰이므로 '
            '(a)에 정의를 그림으로 박아 두었다.',
        without='이 그림이 없으면 논문은 “배아에서 시계가 내려간다”는 기존 관찰의 반복에 그친다. '
                'ZGA 의존성 주장 전체가 이 한 그림에 걸려 있다.',
        panels=[('a', '설계 격자 — 어느 arm에서 언제부터 교란이 걸려 있는지(색칠된 칸), 어느 단계에서 라이브러리를 '
                      '땄는지(점), 그리고 D·I·R의 정의.'),
                ('b', '라이브러리별 시계 값 — 빈 원은 초기 2세포기, 찬 원은 후기 2세포기. arm 위의 숫자가 D, '
                      '교란 arm 아래 숫자가 I와 95% 부트스트랩 구간.')],
        nums=['A485: 대조 −0.24 → −0.06, I = +0.18 [0.11, 0.25]',
              '모계 Brg1 제거 I = +0.11 [−0.01, 0.22]',
              '모계 Tardbp 제거 I = +0.04 (동적 유전자 제거 시 −0.01 → 사전 규칙 미충족)',
              'QC로 76개 중 3개 제외'],
        note='<b>그림을 바꾼 이유:</b> 예전 버전은 막대/forest 형태였다. 그 형태는 (i) 요즘 자동 생성 그림의 전형이라 '
             '인상이 나쁘고, (ii) 그룹당 라이브러리가 2~4개라는 사실을 막대 뒤에 숨긴다. '
             '지금은 점을 전부 보여주어 n을 숨기지 않는다.',
    ),
    dict(
        key='F3', num='그림 3', role='무엇으로 이루어졌나', rolecls='r3',
        q='그 시계 값 하나는 대체 무엇으로 이루어져 있고, ZGA를 되살리면 어떻게 되는가?',
        what='시계 값의 차이를 유전자 1,839개의 기여로 정확히 쪼갠 것. 기여 = (시계 계수) × (그 유전자 발현 변화). '
             '이 항들을 전부 더하면 보고된 시계 차이와 <b>소수점 아래 9자리까지</b> 같다 — 추정이 아니라 항등식이다.',
        why='논문의 심장. 여기서 두 가지 주장이 나온다. ① 보고된 −0.24는 −1.40과 +1.17이라는 큰 두 힘의 '
            '<b>17% 잔차</b>다. ② A485가 없앤 것은 바로 그 잔차를 지배하던 유전자들이다. '
            '그리고 (f)에서 이 논문에서 가장 중요한 장면이 나온다 — DUX로 구제하면 하향 기여는 85% 돌아오는데 '
            '상향 기여도 113%로 같이 커져서, <b>끝점(시계 값)은 제자리에 머문다</b>.',
        without='이 그림이 없으면 “시계 값이 내려갔다/덜 내려갔다”는 현상 기술에서 끝난다. '
                '또한 2.5절의 역설(DUX가 접합자 전사를 되살렸는데 시계는 안 움직임)을 설명할 방법이 사라진다.',
        panels=[('a', '위: 유전자별 기여를 가장 음수부터 가장 양수까지 줄 세운 것. 아래: 그 순서대로 누적합. '
                      '<b>U자가 맞다</b> — 음수부터 더하니 −1.40까지 내려갔다가, 양수 구간에서 올라와 −0.24에서 끝난다. '
                      '그 끝점이 논문이 보고하는 값이다.'),
                ('b', '하향 기여 상위 20개 유전자(검정)와 <b>같은 유전자</b>의 A485에서의 기여(앰버). 오른쪽은 각 유전자의 '
                      'log2 발현 변화. A485가 특정 유전자들의 기여를 지웠다는 것을 보여주는 유일한 패널.'),
                ('c', '두 독립 연구(GSE280522 vs GSE300734)의 유전자별 기여 산점도, r = 0.74. '
                      '“시계 가중치는 잡음”이라는 공격에 대한 방어.'),
                ('d', '접합자 유전자 13개 + 2세포기/DUX 표적 18개의 발현 히트맵. DUX가 실제로 접합자 프로그램을 '
                      '되살렸는지 눈으로 확인. 마커 블록이 대조군 초기부터 이미 높다는 한계도 그림에 그대로 보인다.'),
                ('e', '접합자 유전자 점수 — (d)를 숫자 하나로. A485+DUX − A485 = +1.67.'),
                ('f', '세 arm의 누적 기여 곡선을 <b>대조군 순서 그대로</b> 쌓은 것(같은 유전자를 비교하기 위해). '
                      '점선 = 대조군에서 음의 기여를 낸 666개 유전자 지점. 거기서 −1.40 / −0.76 / −1.19를 읽고, '
                      '곡선의 끝점이 각 arm의 시계 값 −0.24 / −0.06 / −0.06이다.')],
        nums=['대조군: 하향 −1.40, 상향 +1.17 → 잔차 −0.24 (17%)',
              '두 독립 연구 기여 상관 r = 0.74 (상위 50개 중 28개 겹침)',
              '구제 arm: 하향 기여 85% 회복, 상향 기여 113%, 순 값 −0.06',
              '최대 기여 유전자 Klf9 (+4.5), Neto2, Pi4k2a, Gpatch4, Psmb5'],
        note='<b>그림을 합친 이유:</b> 원래 그림 3(분해)과 그림 4(구제)는 따로였다. 둘은 같은 분해 도구를 쓰고 '
             '같은 데이터를 보므로, 합쳐야 “기여가 사라졌다 → 구제하면 돌아온다 → 그래도 순 값은 그대로”가 '
             '한 장 안에서 읽힌다. 본문 그림도 5개에서 4개로 줄어 Aging Cell 제한(6개)에 여유가 생겼다.',
    ),
    dict(
        key='F4', num='그림 4', role='다른 자로 재보기', rolecls='r4',
        q='시계 가중치를 전혀 쓰지 않는 두 번째 자로 재도 같은 답이 나오는가?',
        what='같은 라이브러리에서 splicing(이어맞추기)이 초기→후기 2세포기에 얼마나 진행되는지. '
             '시계 계수가 한 번도 등장하지 않는 독립 지표다.',
        why='시계 지표의 가장 큰 약점은 “가중합 하나”라는 점이다. 그래서 완전히 다른 원리의 자를 하나 더 댄다. '
            '결과적으로 <b>splicing 지표는 사전 규칙을 충족했고 시계 지표는 못 했다</b> — 이 대비 자체가 논문의 '
            '메시지(스칼라 시계는 둔감하다)를 강화한다. (d)는 우리가 스스로 찾아낸 편향을 공개하는 패널이다.',
        without='이 그림이 없으면 모든 결론이 시계 하나에 의존한다. 또 (d)가 없으면 “대조군으로 이벤트를 고르고 '
                '그 대조군을 다시 채점했다”는 순환 논리 공격을 막을 수 없다.',
        panels=[('a', '대조군에서 정의한 splicing 변화 이벤트 511개의 PSI 히트맵. 라이브러리 하나가 열 하나.'),
                ('b', 'arm별 ΔPSI 누적분포(CDF) — splicing 논문의 표준 형태. 중앙값 0.32(대조) / 0.16(A485) / '
                      '0.25(구제). 대조군은 선택 때문에 치우쳐 있다는 점을 그림 설명에 명시했다.'),
                ('c', '진행 점수 — 각 라이브러리를 “대조군 초기 = 0, 대조군 후기 = 1” 축 위에 올린 것. '
                      '세 데이터셋을 한 눈에 비교한다.'),
                ('d', '교차 검증 — 대조군 라이브러리를 하나씩 빼고 다시 정의·채점한 16개 fold. '
                      '표본 내 값의 약 0.21이 선택에서 온 것이었고, 상호작용은 −0.44에서 −0.27로 줄었다. '
                      '16개 중 14개만 음수라 사전 기준에 하나 모자랐고, 그래서 <b>크기는 쓰지 않고 방향만</b> 쓴다.')],
        nums=['대조군 정의 이벤트 511 / 1,376 / 77개 (3개 중 2개에서 기준 충족)',
              '진행 상호작용 −0.44 / −0.42 / −0.45 (세 데이터)',
              '교차 검증 후 −0.27 (16 fold 중 14개 음수)',
              '구제 차이 +0.22 (16 fold 전부 양수)'],
    ),
    dict(
        key='F5', num='그림 5', role='사람에서도 같은가', rolecls='r5',
        q='이 "거의 상쇄"는 생쥐 2세포기만의 성질인가, 시계라는 도구의 성질인가?',
        what='사람 배아 공개 데이터(GSE36552)를 배아 20개로 묶어 같은 시계를 대고, 같은 분해를 한 것. '
             '구간은 결과를 보기 전에 8세포기 → 상실배로 잠갔다(선행연구가 사람에서 하락을 보고한 구간).',
        why='논문 전체가 생쥐 한 종·한 구간에 걸려 있다는 것이 심사에서 가장 먼저 나올 약점이었다. '
            '이 그림은 두 가지를 동시에 한다. ① <b>17% 잔차라는 구조가 다른 종·다른 구간에서도 나온다</b>는 것 '
            '(사람 8%) — 게다가 사람의 그 구간은 <b>세포분열을 건너뛴다</b>. 생쥐 무대의 특징이 없는데도 같은 '
            '모양이 나왔으므로, 이 성질은 2세포기의 특성이 아니라 가중치 덧셈식이 전사체 개편을 읽는 방식에 '
            '가깝다. ② 그런데 <b>유전자 명단은 공유되지 않는다</b>(r = 0.07). 같은 구간을 본 생쥐 연구 두 건이 '
            'r = 0.74였으니 계산 잡음이 아니다.',
        without='이 그림이 없으면 “한 종, 한 구간, 한 시계”라는 비판에 답할 것이 없고, 동시에 기여도 상위 '
                '유전자를 생물학적 주역처럼 읽으려는 유혹을 막을 근거도 없다. 두 주장은 같은 데이터에서 나온다.',
        panels=[('a', '사람 배아 20개의 시계 값(난자 기준). 연한 띠가 사전에 잠근 구간. '
                      '그 구간에서 −0.10으로 궤적 중 가장 큰 하락이지만, 상실배가 2개뿐이라 구간 추정은 0을 '
                      '포함한다 — 그림과 본문에 그대로 적었다.'),
                ('b', '그림 3a와 <b>같은 계산</b>을 사람에서 한 것. 696개가 −1.24까지 끌어내리고, 674개가 '
                      '+1.15로 되밀어, 보고되는 −0.10이 남는다. 잔차 비율 8%.'),
                ('c', '유전자별 기여를 사람(세로) 대 생쥐(가로)로 찍은 산점도. 두 축 모두 기여가 0이 아닌 '
                      '1,098개 유전자. 대각선에 몰리지 않는다 — r = 0.07.')],
        nums=['시계 유전자 1,839개 중 1,370개 검출 (0.74, 생쥐 종간 데이터의 912개보다 많음)',
              '8세포기 → 상실배 −0.098 [−0.23, +0.04]',
              '하향 696개 −1.24 / 상향 674개 +1.15 → 잔차 8% (생쥐 17%)',
              '사람 대 생쥐 기여 상관 r = 0.07 (생쥐 대 생쥐는 0.74)'],
        note='<b>계획서를 고친 기록:</b> 동결 계획서에는 "사람은 시계의 훈련 종이므로 상동유전자 대응이 필요 '
             '없다"고 적혀 있었는데, 실제로는 시계의 입력 10,487개가 전부 <b>생쥐</b> Entrez 번호였다'
             '(사람 표에서는 0개 일치). 첫 실행이 빈 행렬로 실패해 <b>점수가 나오기 전에</b> 이 사실을 '
             'addendum으로 적고, 시계 패키지가 함께 배포하는 대응표로 세 단계 매핑을 선언한 뒤 다시 돌렸다. '
             '판정 규칙과 구간은 바꾸지 않았다.',
    ),
]

SUPP_FIGS = [
    ('S1', '보충 그림 S1', '측정한 ZGA 시점과 시계 (사후)',
     '문헌에서 가져온 major ZGA 구간 대신, <b>같은 배아에서 직접 측정한</b> 접합자 전사체 증가 시점을 썼을 때의 결과. '
     '구간 일치는 여전히 4종 중 1종이지만, 접합자 증가와 시계 하락은 구간에 걸쳐 상관한다(ρ = 0.52). '
     '본문이 아니라 보충인 이유: 사전 계획에 없던 사후 분석이라 본문 주장의 근거로 쓰지 않는다.'),
    ('S2', '보충 그림 S2', '두 번째 bulk 데이터 (GSE66582)',
     '세 번째 재현. 단계당 라이브러리가 2개뿐이라 방향만 읽는다. 크기를 주장할 수 없어서 보충으로 내렸다.'),
    ('S3', '보충 그림 S3', '라이브러리 품질 관리 (QC)',
     '76개 라이브러리 전부의 정렬률과 검출 유전자 수, 그리고 사전 기준선. '
     '<b>제외한 3개가 어느 것인지 이름까지 적혀 있다.</b> "마음에 안 드는 표본을 뺐다"는 의심을 원천 차단하는 그림.'),
    ('S4', '보충 그림 S4', '보조 데이터의 단일 시점 비교',
     'α-amanitin, Obox3 knockdown, 핵 이식(SCNT) 결과. 이 세 가지는 초기/후기 2세포기 쌍이 아니라 한 시점만 '
     '있어서 D를 계산할 수 없다. 그래서 주요 데이터가 아니라 보조로 분리했고, 본문에서는 방향만 인용한다.'),
    ('S5', '보충 그림 S5', '최대 기여 유전자 5개의 발현 + 기여 프로필',
     '(a) Klf9·Neto2·Pi4k2a·Gpatch4·Psmb5가 라이브러리마다 실제로 어떻게 변하는지. 그림 3b의 "기여"가 '
     '추상적 계산값이 아니라 실제 발현 변화에서 나온다는 뒷받침. (b) 유전자별 기여를 교란·구제 arm 대 대조군으로 '
     '찍은 산점도 — 대조군과의 상관이 A485에서 0.65, 구제 arm에서 0.83으로 올라간다. '
     '<b>원래 보충 그림 두 장이었는데 같은 데이터의 사후 상세라 하나로 합쳤다.</b>'),
    ('S6', '보충 그림 S6', 'splicing 이벤트 선택의 효과크기와 p 값 (volcano)',
     '필터를 통과한 이벤트 전부를 가로축 ΔPSI(대조군 후기 − 초기), 세로축 −log10 p로 찍고, 사전에 정한 선택 '
     '기준(|ΔPSI| ≥ 0.10, 유전자 보정 empirical p &lt; 0.05)을 점선으로 그렸다. 검은 점이 선택된 511 / 1,376 / 77개. '
     '이 p 값은 선택에만 쓰고 추론에는 쓰지 않는다.'),
    ('S7', '보충 그림 S7', 'junction read로 splicing을 다시 재본 검증 (사후)',
     '<b>새 분석.</b> GSE280522 23개 라이브러리를 STAR로 다시 정렬해 junction read를 직접 세었다. '
     '(a) junction PSI와 전사체 PSI의 밀도 지도 — 58,958쌍에서 r = ρ = 0.78로 일치한다. '
     '(b) 규칙으로 뽑은 3개 이벤트의 read 커버리지와 junction 수: A485에서 inclusion junction이 떨어지고 '
     'DUX 구제에서 돌아오는 것이 read 수준에서 그대로 보인다. '
     '논문의 한계 3번("junction read가 아니라 전사체 추정치로 정량했다")을 직접 검정한 그림이다.'),
]

GATES = [
    ('1', '가장 큰 하락이 문헌 major ZGA 구간에 있는가?', '4종 중 3종 이상', '4종 중 1종 — <b>미충족</b>', 'bad'),
    ('1-alt', '(사후) 직접 측정한 ZGA 시점을 쓰면?', '구간 일치 + 상관', '구간은 1/4, 상관은 있음(ρ 0.52)', 'mid'),
    ('2a', '모계 RNA 청소만으로 생기는 착시인가?', '잔존율과 모의 재현율', '약 절반은 청소로 재현 — <b>중간</b>', 'mid'),
    ('2b', '2세포기 안 하락이 다른 데이터에서도 보이는가?', '후기 − 초기 < 0', '−0.14, −0.26 — <b>재현</b>', 'good'),
    ('3', 'ZGA를 막으면 하락이 약해지는가?', '판정 가능한 주요 데이터 전부에서 I > 0, 두 유전자 집합 일치',
     '전체 유전자에서는 성립, 동적 유전자를 빼면 가장 작은 데이터 하나가 뒤집힘 — <b>완전 충족은 아님</b>', 'mid'),
    ('4', 'splicing 진행도 ZGA에 의존하는가?', '3개 중 2개에서 음의 상호작용 + 양의 구제',
     '<b>충족</b> (크기는 교차 검증으로 제한)', 'good'),
    ('H', '(사후, 점수 전 동결) 사람 배아도 같은 구조인가?', '잔차 비율이 생쥐와 같은 자릿수',
     '8% 대 17% — <b>같은 구조</b> (단, 유전자는 다름)', 'good'),
]

REVISED = [
    ('논문의 중심을 옮긴 것', [
        '제목과 초록의 중심을 <b>ZGA 의존성</b>에서 <b>“하락은 거의 상쇄되고 남은 나머지”</b>로 옮겼다. '
        '이유는 근거의 세기다 — ZGA 의존은 사실상 데이터 한 건(GSE280522)이 떠받치고 있는데, 분해와 은폐는 '
        '데이터 세 건과 두 종에서 같은 모양으로 나온다. 약한 쪽을 제목에 걸면 심사에서 거기부터 무너진다.',
        '그래서 <b>사람 배아 재현(그림 5)</b>을 새로 붙였다. 다른 종, 다른 구간, 분열을 건너뛰는 구간에서도 '
        '같은 구조가 나오는지 보기 위해서다 — 나왔다(8% 대 17%).',
        'ZGA 교란(그림 2·3·4)은 주장에서 <b>무대</b>로 내렸다. 결과는 그대로 두고, 그것이 무엇을 보이는지의 '
        '위치만 바꿨다.',
        '동시에 알게 된 것 하나 — 기여도 상위 유전자는 종을 넘으면 남지 않는다(r = 0.07). 그래서 그 목록을 '
        '생물학적 발견처럼 읽지 말라는 문장을 고찰에 넣었다.',
    ]),
    ('말이 과했던 것', [
        '제목과 소제목의 “minor ZGA” → “ZGA”. A485는 minor와 major를 구분해서 막지 못하므로 데이터가 '
        '“minor”를 집어낼 수 없다. 대신 왜 구분이 안 되는지 한 문장을 본문에 넣었다.',
        '“largely restores(대부분 회복시킨다)” 같은 표현을 실제 수치(하향 85%, 접합자 전사 66%)로 바꿨다.',
        '“splicing이 더 깨끗한 답을 줬다”는 문장에 단서를 달았다 — 두 지표는 검정력도 데이터도 다르다.',
        '“하락의 약 절반이 모계 청소” → 0.70(생쥐), 0.47(소)로 숫자를 그대로 적었다.',
    ]),
    ('숫자의 기준을 명시', [
        '85% / 98% / 113%가 각각 무엇에 대한 비율인지 본문에 적었다(모두 대조군 대비).',
        '교차 검증 fold 범위를 신뢰구간이 아니라 “흩어짐”으로 부른다고 방법에 명시했다.',
        '상호작용 +0.18 자체가 이미 잔차(+0.64 − 0.46)라는 점을 밝혔다.',
        'Klf9 등 개별 유전자가 A485에서 기여의 몇 %를 잃었는지 적었다(47~88%).',
    ]),
    ('빠졌던 것을 채움', [
        '그림 3d 히트맵에 Zfp352가 들어 있는데 설명에는 빠져 있었다 → 설명과 방법의 마커 목록에 추가.',
        '데이터별 read 길이, QC 규칙, 난수 seed, SUPPA2 설정을 본문/보충에 명시.',
        'GSE225056의 5′ 행렬이 361개 중 340개만 담고 있다는 사실을 적었다.',
    ]),
    ('그림 쪽', [
        '보충 그림 번호를 <b>본문에서 처음 인용되는 순서</b>로 다시 매겼다(예전에는 순서가 뒤섞여 있었다).',
        '모든 그림의 글자 크기를 6 pt 이상으로 올렸다(저널 최소 기준).',
        '그림 2의 각주와 라벨을 실제 본문 위치와 맞췄다.',
        '<b>그림 11개를 전수 검사해 캔버스 밖으로 잘린 글자를 6곳에서 찾아 고쳤다</b> — 그림 1d·2a의 각주, '
        '그림 3f·4b의 축 라벨, 보충 S1·S4의 각주. 인쇄하면 문장이 잘린 채 나가는 상태였다.',
    ]),
]

FACTCHECK = [
    ('원문과 대조해 고친 것', [
        '<b>DBTMEE 유전자 분류.</b> "다섯 개 모두 major ZGA·2세포기 일과성·전환기"라고 썼는데, 데이터베이스 원본 표를 '
        '직접 내려받아 조회하니 <b>Psmb5는 MGA</b>(착상 전 중기 물결)였고 "전환기"는 DBTMEE의 분류명이 아니었다. '
        '이제 유전자별로 실제 분류명을 적는다.',
        '<b>존재하지 않는 인용.</b> Tyshkovskiy 2026에서 "MYC 표적 프로그램"을 인용했으나 그 논문에 없는 내용이었고, '
        '세포주기 모듈은 오히려 초기 발생에서 tAge가 <b>올라간다</b>고 적혀 있었다 — 우리가 쓴 것과 반대.',
        '<b>돼지 음성 주장.</b> "돼지 착상 전 배아에 분자 나이 지표가 없다"고 썼으나 텔로미어 측정이 존재했다'
        '(Dang-Nguyen 2012). 텔로미어 신장을 "두 종"이라 한 것도 틀렸다.',
        '<b>리뷰를 1차 문헌처럼 인용.</b> 텔로미어 신장의 출처를 2005년 리뷰로 달았으나 1차 문헌은 '
        'Schaetzlein et al. 2004 <i>PNAS</i>다.',
        '<b>인용 대상 착오.</b> Higgins-Chen 2022는 시계 <b>가중치</b>가 아니라 개별 <b>CpG 측정치</b>의 잡음을 보인 '
        '논문이다(반복 측정 간 최대 9년 편차).',
        '<b>데이터 설계 오기.</b> GSE162345의 핵 융합 시점을 "21–24 / 30–33 hpi"로 썼으나 원문은 <b>21 hpi 또는 '
        '30 hpi 단일 시점</b>이다. 또 그 2세포기 arm이 원 연구에서는 4세포기 시스템의 <b>대조군</b>이고 30 hpi에서는 '
        '공여 유전체 유전자가 8개만 활성화됐다는 사실을 본문에 넣었다.',
        '<b>데이터 성격 오기.</b> GSE248499를 "히스톤 탈메틸화효소 핵 이식"이라 불렀으나 원 논문은 <b>G9a 저해제</b> '
        '연구다(G9a는 메틸전이효소). 우리가 실제로 쓴 arm만 기술하도록 고쳤다.',
        '<b>엉뚱한 인용.</b> "종간 ZGA 비교"의 근거로 든 Li et al. 2025는 SCNT 재프로그래밍 비교였다. 정답은 우리가 '
        '이미 쓰고 있던 Oomen et al. 2025였다.',
        '<b>지지하지 않는 인용.</b> Hou et al. 2025는 전문에 "상피화"가 한 번도 나오지 않는다. 해당 문장은 '
        'Chandramohan 2026만 남겼다.',
    ]),
    ('논증의 강도를 낮춘 것', [
        '<b>A485는 전사 저해제가 아니다.</b> Xiao et al. 2025에서 A485의 ZGA 효과는 <i>Dux</i> 유도 실패를 거쳐 '
        '나타나고 외인성 Dux로 우회된다. 따라서 "아세틸화 의존 ZGA"를 보는 지표이지 <b>새로운 전사에 대한 의존의 '
        '독립 확증이 아니다.</b> 1차 약리 문헌(Lasko 2017)과 배아 약리 문헌(Wang 2022)을 함께 인용했다.',
        '<b>Obox3 다리도 약하다.</b> 원 논문의 Obox3 근거는 전부 획득형(배아줄기세포 과발현, SCNT 회복)이고 '
        'knockdown 결과는 초록에 없다. OBOX의 ZGA 요구성도 6중 녹아웃에서 보인 것이다.',
        '결과적으로 "ZGA 의존"의 세 다리 중 <b>Pol II를 직접 저해하는 것은 α-amanitin뿐</b>임을 고찰에 명시했다.',
    ]),
    ('종간 비교에 추가한 한계', [
        '생쥐 구간만 <b>분열이 없는</b> 구간이고 사전에 정한 돼지·소·토끼 구간은 모두 분열을 하나씩 포함한다 — '
        '네 검정은 같은 검정이 아니었다.',
        '돼지는 major 활성화가 <b>4세포기 안</b>에서 일어나 생쥐와 같은 구조를 만들 수 있었지만'
        '(Prather & Rickords 1992), 등록 자료에 채취 시점이 없는 4세포기 시료 한 묶음뿐이라 단계 내부 대조를 '
        '만들 수 없었다.',
        '따라서 "4종 중 1종"이라는 음성 결과는 시계의 실패가 아니라 <b>설계·데이터의 한계</b>로도 설명된다.',
    ]),
]

KEPT = [
    '그룹당 라이브러리가 2~4개다. 더 큰 공개 데이터가 없어 해결이 불가능하며, 그래서 p 값 대신 효과 크기·구간·'
    '데이터 간 방향 일치로 논증한다. (한계 첫 번째로 이미 적혀 있다)',
    '교란마다 ZGA에 작용하는 경로가 다르고 ZGA 이외의 효과도 있다. 약물·유전자 제거가 ZGA만 건드린다고 '
    '주장하지 않는다.',
    '핵 이식(SCNT) 배아가 체외수정보다 낮게 채점되는 이유는 설명하지 못한다. <b>사전 예상과 반대</b>이고, '
    '그대로 적은 뒤 관련 문헌(Matoba 2014/2024, 복제 동물 수명 보고)을 인용해 부호를 열어 두었다.',
    '교차 검증은 대조군 라이브러리가 충분한 한 데이터에서만 가능했다.',
    'bulk splicing 점수는 조성 변화와 세포 내부 변화를 혼동할 수 있다(한 단계의 배아 전체를 쓰므로 줄지만 '
    '사라지지는 않는다).',
    '분해·구제·교차 검증은 모두 사후 분석이다. 각각 실행 전에 계획을 동결했지만 gate 판정을 본 뒤였다는 점을 '
    '모든 해당 절 제목에 “(사후)”로 표시했다.',
]

# ------------------------------------------------------------------ html ----

def panel_list(panels):
    return '\n'.join(
        f'<li><span class="pk">{k}</span><span class="pt">{t}</span></li>' for k, t in panels)


def fig_card(f):
    note = f'<p class="note">{f["note"]}</p>' if f.get('note') else ''
    nums = '\n'.join(f'<li>{n}</li>' for n in f['nums'])
    return f'''
<section class="figcard {f['rolecls']}" id="fig{f['key']}">
  <div class="fighead">
    <span class="eyebrow">{f['num']}</span>
    <span class="chip">{f['role']}</span>
  </div>
  <h3>{f['q']}</h3>
  <figure><img src="{img(f['key'])}" alt="{f['num']}"></figure>
  <div class="qa">
    <div class="qa-row"><div class="qa-k">무엇을 보나</div><div class="qa-v">{f['what']}</div></div>
    <div class="qa-row hi"><div class="qa-k">왜 필요한가</div><div class="qa-v">{f['why']}</div></div>
    <div class="qa-row"><div class="qa-k">없으면</div><div class="qa-v">{f['without']}</div></div>
  </div>
  <div class="panels">
    <h4>패널별로</h4>
    <ul>{panel_list(f['panels'])}</ul>
  </div>
  <div class="nums"><h4>핵심 숫자</h4><ul>{nums}</ul></div>
  {note}
</section>'''


def supp_card(k, num, title, text):
    return f'''
<section class="suppcard">
  <div class="supptext">
    <span class="eyebrow">{num}</span>
    <h4>{title}</h4>
    <p>{text}</p>
  </div>
  <figure><img src="{img(k)}" alt="{num}"></figure>
</section>'''


gate_rows = '\n'.join(
    f'<tr><td class="g-id">{g[0]}</td><td>{g[1]}</td><td class="g-rule">{g[2]}</td>'
    f'<td class="g-{g[4]}">{g[3]}</td></tr>' for g in GATES)

fact_blocks = '\n'.join(
    f'<div class="revblock"><h4>{ti}</h4><ul>' + '\n'.join(f'<li>{i}</li>' for i in items) + '</ul></div>'
    for ti, items in FACTCHECK)

rev_blocks = '\n'.join(
    f'<div class="revblock"><h4>{t}</h4><ul>' + '\n'.join(f'<li>{i}</li>' for i in items) + '</ul></div>'
    for t, items in REVISED)

kept_items = '\n'.join(f'<li>{k}</li>' for k in KEPT)

HTML = f'''<title>배아 시계 하락, 쉽게 읽기</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Gowun+Batang:wght@400;700&family=Noto+Sans+KR:wght@300;400;500;700&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
:root {{
  --paper:#FAF8F4; --card:#FFFDF9; --ink:#211E1A; --ink2:#413C34; --muted:#6F6A61;
  --rule:#E4DFD5; --rule2:#D6D0C3;
  --amber:#A8620A; --teal:#0F6F68; --down:#9A3A2E; --up:#2F5B86; --violet:#5C4B8A;
  --good:#2F6B43; --mid:#8A6A12; --bad:#8E3A2C;
  --shadow:0 1px 2px rgba(33,30,26,.05), 0 8px 24px -16px rgba(33,30,26,.25);
}}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{
    --paper:#15130F; --card:#1C1A15; --ink:#EEEAE1; --ink2:#D4CEC2; --muted:#9A9488;
    --rule:#2E2A23; --rule2:#3C372E;
    --amber:#D89440; --teal:#4FB3A8; --down:#D2766A; --up:#7FA8CF; --violet:#A192CE;
    --good:#72B183; --mid:#C8A34C; --bad:#D2837A;
    --shadow:0 1px 2px rgba(0,0,0,.4), 0 8px 24px -16px rgba(0,0,0,.8);
  }}
}}
:root[data-theme="dark"] {{
  --paper:#15130F; --card:#1C1A15; --ink:#EEEAE1; --ink2:#D4CEC2; --muted:#9A9488;
  --rule:#2E2A23; --rule2:#3C372E;
  --amber:#D89440; --teal:#4FB3A8; --down:#D2766A; --up:#7FA8CF; --violet:#A192CE;
  --good:#72B183; --mid:#C8A34C; --bad:#D2837A;
  --shadow:0 1px 2px rgba(0,0,0,.4), 0 8px 24px -16px rgba(0,0,0,.8);
}}
* {{ box-sizing:border-box; }}
body {{
  margin:0; background:var(--paper); color:var(--ink);
  font-family:'Noto Sans KR', -apple-system, 'Apple SD Gothic Neo', 'Malgun Gothic', sans-serif;
  font-weight:400; font-size:15.5px; line-height:1.85; word-break:keep-all;
  -webkit-font-smoothing:antialiased;
}}
.wrap {{ max-width:1000px; margin:0 auto; padding-block:0 96px; padding-left:20px; padding-right:20px; }}
.measure {{ max-width:760px; }}
h1,h2,h3,h4 {{ font-family:'Gowun Batang', 'Noto Serif KR', Georgia, serif; text-wrap:balance; margin:0; }}
p {{ margin:0 0 .9em; }}
b {{ font-weight:700; color:var(--ink); }}
a {{ color:inherit; }}

/* header */
header.top {{
  position:sticky; top:0; z-index:20; background:color-mix(in srgb, var(--paper) 92%, transparent);
  backdrop-filter:blur(8px); border-bottom:1px solid var(--rule);
}}
.topin {{ max-width:1000px; margin:0 auto; padding:10px 20px; display:flex; gap:14px;
  align-items:baseline; flex-wrap:wrap; }}
.topin .brand {{ font-family:'Gowun Batang',serif; font-weight:700; font-size:14px; white-space:nowrap; }}
nav {{ display:flex; gap:4px; flex-wrap:wrap; }}
nav a {{ font-size:12px; color:var(--muted); text-decoration:none; padding:3px 8px; border-radius:999px;
  border:1px solid transparent; }}
nav a:hover {{ color:var(--ink); border-color:var(--rule2); }}

/* hero */
.hero {{ padding-block:56px 40px; border-bottom:1px solid var(--rule); }}
.hero .kicker {{ font-family:'IBM Plex Mono',monospace; font-size:11.5px; letter-spacing:.12em;
  text-transform:uppercase; color:var(--amber); margin-bottom:18px; }}
.hero h1 {{ font-size:clamp(28px,4.6vw,42px); line-height:1.35; margin-bottom:22px; }}
.hero .lede {{ font-size:17px; color:var(--ink2); max-width:720px; }}
.oneline {{ margin-top:28px; padding:18px 22px; border-left:3px solid var(--amber);
  background:var(--card); border-radius:0 10px 10px 0; box-shadow:var(--shadow); }}
.oneline .l {{ font-family:'IBM Plex Mono',monospace; font-size:11px; letter-spacing:.1em;
  color:var(--muted); text-transform:uppercase; display:block; margin-bottom:6px; }}
.oneline p {{ margin:0; font-size:16.5px; }}

h2.sec {{ font-size:25px; margin:64px 0 6px; }}
h2.sec + .sub {{ color:var(--muted); font-size:14px; margin-bottom:26px; }}

/* problem cards */
.probs {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(230px,1fr)); gap:16px; }}
.prob {{ background:var(--card); border:1px solid var(--rule); border-radius:12px; padding:20px 20px 18px; }}
.prob .n {{ font-family:'IBM Plex Mono',monospace; font-size:11px; color:var(--amber); letter-spacing:.1em; }}
.prob h4 {{ font-size:16.5px; margin:8px 0 8px; }}
.prob p {{ font-size:14px; color:var(--ink2); margin:0; line-height:1.75; }}

/* chain */
.chain {{ counter-reset:step; display:flex; flex-direction:column; gap:0; }}
.step {{ position:relative; padding:18px 0 18px 58px; border-bottom:1px dashed var(--rule); }}
.step:last-child {{ border-bottom:none; }}
.step::before {{
  counter-increment:step; content:counter(step);
  position:absolute; left:0; top:20px; width:32px; height:32px; border-radius:50%;
  display:grid; place-items:center; font-family:'IBM Plex Mono',monospace; font-size:13px;
  background:var(--card); border:1px solid var(--rule2); color:var(--amber);
}}
.step h4 {{ font-size:17px; margin-bottom:5px; }}
.step p {{ font-size:14.5px; color:var(--ink2); margin:0 0 6px; }}
.step .where {{ font-family:'IBM Plex Mono',monospace; font-size:11.5px; color:var(--muted); }}

/* figure cards */
.figcard {{
  background:var(--card); border:1px solid var(--rule); border-radius:14px;
  padding:26px 26px 22px; margin:26px 0; box-shadow:var(--shadow);
  border-top:3px solid var(--rule2);
}}
.figcard.r1 {{ border-top-color:var(--up); }}
.figcard.r2 {{ border-top-color:var(--amber); }}
.figcard.r3 {{ border-top-color:var(--down); }}
.figcard.r4 {{ border-top-color:var(--teal); }}
.figcard.r5 {{ border-top-color:var(--violet); }}
.fighead {{ display:flex; align-items:center; gap:10px; margin-bottom:10px; flex-wrap:wrap; }}
.eyebrow {{ font-family:'IBM Plex Mono',monospace; font-size:11.5px; letter-spacing:.12em;
  text-transform:uppercase; color:var(--muted); }}
.chip {{ font-size:11.5px; padding:2px 10px; border-radius:999px; border:1px solid var(--rule2);
  color:var(--ink2); background:var(--paper); }}
.figcard h3 {{ font-size:21px; line-height:1.5; margin-bottom:18px; }}
.figcard figure {{ margin:0 0 20px; background:#fff; border:1px solid var(--rule); border-radius:8px;
  padding:10px; overflow-x:auto; }}
.figcard figure img {{ display:block; width:100%; max-width:100%; height:auto; }}

.qa {{ display:flex; flex-direction:column; gap:2px; margin-bottom:20px; }}
.qa-row {{ display:grid; grid-template-columns:86px 1fr; gap:16px; padding:11px 0;
  border-bottom:1px solid var(--rule); }}
.qa-row:last-child {{ border-bottom:none; }}
.qa-row.hi {{ background:linear-gradient(90deg, color-mix(in srgb, var(--amber) 7%, transparent), transparent 70%);
  border-radius:6px; padding-left:12px; margin-left:-12px; }}
.qa-k {{ font-family:'IBM Plex Mono',monospace; font-size:11.5px; letter-spacing:.06em;
  color:var(--muted); padding-top:5px; }}
.qa-v {{ font-size:14.8px; color:var(--ink2); }}

.panels, .nums {{ margin-top:18px; }}
.panels h4, .nums h4 {{ font-size:12px; font-family:'IBM Plex Mono',monospace; letter-spacing:.1em;
  text-transform:uppercase; color:var(--muted); margin-bottom:10px; font-weight:500; }}
.panels ul {{ list-style:none; margin:0; padding:0; display:flex; flex-direction:column; gap:9px; }}
.panels li {{ display:grid; grid-template-columns:22px 1fr; gap:10px; font-size:14.3px; color:var(--ink2); }}
.pk {{ font-family:'IBM Plex Mono',monospace; font-weight:500; color:var(--ink);
  background:var(--paper); border:1px solid var(--rule2); border-radius:5px; height:21px;
  display:grid; place-items:center; font-size:11.5px; margin-top:4px; }}
.nums ul {{ list-style:none; margin:0; padding:0; display:flex; flex-wrap:wrap; gap:7px; }}
.nums li {{ font-family:'IBM Plex Mono',monospace; font-size:12.2px; color:var(--ink2);
  background:var(--paper); border:1px solid var(--rule); border-radius:6px; padding:5px 10px; }}
.note {{ margin-top:18px; padding:14px 16px; border-left:2px solid var(--violet);
  background:color-mix(in srgb, var(--violet) 6%, transparent); border-radius:0 8px 8px 0;
  font-size:14px; color:var(--ink2); }}

/* supp */
.suppcard {{ display:grid; grid-template-columns:1fr 340px; gap:22px; align-items:start;
  padding:22px 0; border-bottom:1px solid var(--rule); }}
.suppcard:last-child {{ border-bottom:none; }}
.suppcard h4 {{ font-size:16.5px; margin:6px 0 8px; }}
.suppcard p {{ font-size:14px; color:var(--ink2); margin:0; }}
.suppcard figure {{ margin:0; background:#fff; border:1px solid var(--rule); border-radius:6px; padding:7px; }}
.suppcard figure img {{ display:block; width:100%; height:auto; }}

/* table */
.tablewrap {{ overflow-x:auto; border:1px solid var(--rule); border-radius:10px; background:var(--card); }}
table {{ border-collapse:collapse; width:100%; min-width:640px; font-size:14px; }}
th, td {{ text-align:left; padding:12px 14px; border-bottom:1px solid var(--rule); vertical-align:top; }}
th {{ font-family:'IBM Plex Mono',monospace; font-size:11px; letter-spacing:.08em; text-transform:uppercase;
  color:var(--muted); font-weight:500; }}
tr:last-child td {{ border-bottom:none; }}
.g-id {{ font-family:'IBM Plex Mono',monospace; color:var(--amber); white-space:nowrap; }}
.g-rule {{ color:var(--muted); font-size:13px; }}
.g-good {{ color:var(--good); }} .g-mid {{ color:var(--mid); }} .g-bad {{ color:var(--bad); }}

/* revision */
.revgrid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(290px,1fr)); gap:18px; }}
.revblock {{ background:var(--card); border:1px solid var(--rule); border-radius:12px; padding:20px; }}
.revblock h4 {{ font-size:16px; margin-bottom:12px; }}
.revblock ul {{ margin:0; padding-left:18px; display:flex; flex-direction:column; gap:9px; }}
.revblock li {{ font-size:14px; color:var(--ink2); }}
.kept {{ background:var(--card); border:1px solid var(--rule); border-left:3px solid var(--bad);
  border-radius:0 12px 12px 0; padding:22px 24px; }}
.kept ul {{ margin:0; padding-left:18px; display:flex; flex-direction:column; gap:11px; }}
.kept li {{ font-size:14.3px; color:var(--ink2); }}

/* files */
.files {{ font-family:'IBM Plex Mono',monospace; font-size:12.5px; }}
.files td:first-child {{ color:var(--muted); }}

footer {{ margin-top:72px; padding-top:22px; border-top:1px solid var(--rule); color:var(--muted);
  font-size:12.5px; }}

@media (max-width:760px) {{
  .suppcard {{ grid-template-columns:1fr; }}
  .qa-row {{ grid-template-columns:1fr; gap:3px; }}
  .qa-k {{ padding-top:0; }}
  .figcard {{ padding:20px 16px 18px; }}
}}
</style>

<header class="top">
  <div class="topin">
    <span class="brand">2세포기 시계 하락 · 쉽게 읽기</span>
    <nav>
      <a href="#why">왜 이 연구</a>
      <a href="#chain">논리 사슬</a>
      <a href="#figs">본문 그림</a>
      <a href="#supp">보충 그림</a>
      <a href="#gates">판정</a>
      <a href="#review">리뷰와 수정</a>
      <a href="#factcheck">문헌 대조</a>
      <a href="#files">파일</a>
    </nav>
  </div>
</header>

<div class="wrap">

<div class="hero measure">
  <div class="kicker">원고 v1 · 쉬운 말 요약 · 2026-09-29 갱신</div>
  <h1>배아의 전사체 나이 하락은 무엇으로 이루어져 있고,<br>각 그림은 왜 거기 있는가</h1>
  <p class="lede">
    배아는 발생 초기에 분자 나이가 <b>젊어진다</b>고 알려져 있다. 이 원고는 그 하락을 생쥐 2세포기 안에서 붙잡아
    시계 값 하나를 <b>유전자 단위로 완전히 분해</b>하고, <b>ZGA를 막았을 때·되살렸을 때</b> 무엇이 움직이는지 보고,
    같은 계산을 <b>사람 배아</b>에서 되풀이한다. 아래는 논문을 처음 보는 사람 기준으로 다시 쓴 설명이다.
  </p>
  <div class="oneline">
    <span class="l">한 문장</span>
    <p>보고되는 −0.24라는 하락은 −1.40과 +1.17이라는 큰 두 힘이 남긴 <b>17%의 잔차</b>다.
       사람 배아의 다른 구간에서도 같은 모양(8%)이 나오지만 <b>유전자는 공유되지 않는다</b>.
       그리고 ZGA를 되살리면 그 기여들은 돌아오는데도 <b>시계 값은 제자리에 머문다</b>.</p>
  </div>
</div>

<h2 class="sec" id="why">왜 이 연구가 필요한가</h2>
<p class="sub">기존 배아 시계 연구가 답하지 못한 세 가지</p>
<div class="probs">
  <div class="prob">
    <div class="n">문제 1</div>
    <h4>관찰만 했다</h4>
    <p>지금까지의 배아 시계 연구는 정상 배아를 따라가며 “내려간다”를 기술했다.
       그 하락이 <b>무엇에 의존하는지</b> 교란 실험으로 검정한 연구는 없었다.</p>
  </div>
  <div class="prob">
    <div class="n">문제 2</div>
    <h4>조성으로도 움직인다</h4>
    <p>시계 값은 세포의 상태가 아니라 <b>표본의 구성</b>이 바뀌어도 움직인다.
       초기 배아에서는 모계 RNA가 대량으로 사라지므로 이 함정이 특히 크다.</p>
  </div>
  <div class="prob">
    <div class="n">문제 3</div>
    <h4>가중합은 상쇄된다</h4>
    <p>시계는 수천 개 유전자의 가중합이다. 구성 요소들이 서로 반대로 크게 움직여도
       <b>합은 거의 움직이지 않을 수 있다</b>. 숫자 하나만 보면 이걸 절대 못 본다.</p>
  </div>
</div>

<h2 class="sec" id="chain">논리 사슬 일곱 단계</h2>
<p class="sub">각 단계가 어느 그림으로 뒷받침되는지</p>
<div class="chain measure">
  <div class="step">
    <h4>하락이 실제로 있는가</h4>
    <p>세포분열이 없는 2세포기 안에서, 독립된 생쥐 데이터 세 건 모두 시계 값이 내려간다.
       종간 비교에서는 사전 기준이 4종 중 1종에서만 충족되었다 — 이 음성 결과도 그대로 보고한다.</p>
    <div class="where">그림 1b, 1c · 보충 S1, S2</div>
  </div>
  <div class="step">
    <h4>단순한 “청소 효과”인가</h4>
    <p>모계 RNA 제거만 모의해도 원래 하락의 0.70(생쥐)·0.47(소)이 재현된다. 절반쯤은 청소 효과다.
       그러나 동적 유전자를 모두 빼도 하락의 일부는 남는다.</p>
    <div class="where">그림 1d</div>
  </div>
  <div class="step">
    <h4>남는 부분이 ZGA에 의존하는가</h4>
    <p>A485로 ZGA를 막으면 하락이 −0.24에서 −0.06으로 약해진다(I = +0.18).
       모계 Brg1 제거·Obox3 knockdown·α-amanitin도 같은 방향. 가장 작은 데이터 하나가 부호를 뒤집어
       사전 규칙은 완전히 충족되지 않았다.</p>
    <div class="where">그림 2b · 보충 S3, S4</div>
  </div>
  <div class="step">
    <h4>그 시계 값은 무엇으로 이루어졌는가</h4>
    <p>유전자별로 정확히 쪼개면(항등식) 하락은 −1.40과 +1.17의 잔차다.
       구성은 독립된 두 연구에서 일치하고(r = 0.74), A485가 없앤 것은 바로 그 지배적 기여들이다.</p>
    <div class="where">그림 3a, 3b, 3c · 보충 S5</div>
  </div>
  <div class="step">
    <h4>되살리면 돌아오는가</h4>
    <p>DUX를 함께 발현시키면 접합자 전사와 하향 기여가 대부분 돌아온다(하향 85%).
       그런데 상향 기여도 113%로 커져서 <b>순 값은 A485 수준에 머문다</b>.
       이것이 “스칼라 하나는 회복을 가린다”의 증거다.</p>
    <div class="where">그림 3d, 3e, 3f · 보충 S6</div>
  </div>
  <div class="step">
    <h4>사람에서도 같은 구조인가</h4>
    <p>사람 배아의 다른 구간(8세포기→상실배, 분열을 건너뛴다)에서도 하락은 −1.24와 +1.15의 <b>8% 잔차</b>다.
       구조는 반복되는데 <b>유전자는 공유되지 않는다</b>(r = 0.07; 생쥐끼리는 0.74).</p>
    <div class="where">그림 5</div>
  </div>
  <div class="step">
    <h4>시계를 안 쓰는 자로도 같은가</h4>
    <p>splicing 진행이라는 독립 지표에서도 ZGA를 막으면 진행이 줄고, 구제하면 대조군 쪽으로 돌아온다.
       이 지표는 사전 규칙을 충족했다 — 시계 지표는 못 했다.</p>
    <div class="where">그림 4 전체 · 보충 S7</div>
  </div>
</div>

<h2 class="sec" id="figs">본문 그림 5개 — 각각 왜 존재하는가</h2>
<p class="sub">그림 하나가 논리 사슬의 한 마디를 맡는다. 하나라도 빠지면 그 자리에서 공격받는다.</p>
{''.join(fig_card(f) for f in MAIN_FIGS)}

<h2 class="sec" id="supp">보충 그림 7개 — 왜 본문이 아니라 보충인가</h2>
<p class="sub">주장을 지탱하지만 본문 논리 사슬의 마디는 아닌 것들. 번호는 본문에서 처음 인용되는 순서다.</p>
{''.join(supp_card(*s) for s in SUPP_FIGS)}

<h2 class="sec" id="gates">사전에 정한 판정과 실제 결과</h2>
<p class="sub">각 분석 전에 규칙을 파일로 동결하고, 데이터를 본 뒤에는 규칙을 바꾸지 않았다.
   그래서 <b>미충족도 그대로 보고한다</b>.</p>
<div class="tablewrap">
  <table>
    <thead><tr><th>Gate</th><th>물음</th><th>사전 규칙</th><th>결과</th></tr></thead>
    <tbody>{gate_rows}</tbody>
  </table>
</div>
<p class="sub" style="margin-top:14px">
  두 개가 완전히 충족되지 않았다는 사실이 이 논문의 약점이자 동시에 신뢰의 근거다.
  규칙을 나중에 고쳤다면 전부 “충족”으로 만들 수 있었지만, 그러지 않았다는 것이 감사 기록에 남아 있다.
</p>

<h2 class="sec" id="review">스스로 한 비판적 리뷰 — 무엇을 지적했고 무엇을 고쳤나</h2>
<p class="sub">리뷰어 관점으로 논리 모순·그림 구성·잘못된 언급을 훑어 66건을 찾았고, 그중 64건이 검증을 통과했다.
   고칠 수 있는 것은 고쳤고, 고칠 수 없는 것은 한계로 남겼다.</p>
<div class="revgrid">{rev_blocks}</div>

<h2 class="sec" id="factcheck">문헌 전수 대조 — 배경 주장을 원문과 맞춰봤다</h2>
<p class="sub">서론과 고찰의 문헌 기반 진술을 인용 논문 원문과 하나씩 대조했다. 실제 오류 20건이 나왔고 전부 고쳤다.
   가장 심각한 두 건(DBTMEE 분류, A485 논리)은 데이터베이스 원본 파일과 1차 논문으로 직접 재확인했다.</p>
<div class="revgrid">{fact_blocks}</div>

<h4 style="font-size:17px; margin:34px 0 12px">고치지 않고 “한계”로 남긴 것</h4>
<div class="kept">
  <ul>{kept_items}</ul>
</div>

<h2 class="sec" id="files">파일은 어디 있나</h2>
<p class="sub">기준 경로: <code>…/embryogenesis_Aging/analysis/</code></p>
<div class="tablewrap">
  <table class="files">
    <tbody>
      <tr><td>원고 (영문·작업본)</td><td>gz/results/MANUSCRIPT_gz_v1_EN.md</td></tr>
      <tr><td>원고 (국문)</td><td>gz/results/MANUSCRIPT_gz_v1_KR.md</td></tr>
      <tr><td>투고 패키지</td><td>gz/results/submission/ (docx · Figure1–4.pdf · FigureS1–S7.pdf · 표)</td></tr>
      <tr><td>본문 그림</td><td>gz/figures/pub/</td></tr>
      <tr><td>보충 그림</td><td>gz/figures/supp/</td></tr>
      <tr><td>그림 원자료</td><td>gz/figures/source_data/</td></tr>
      <tr><td>보충표</td><td>gz/results/supp_tables/</td></tr>
      <tr><td>이 페이지</td><td>gz/results/review/easy_review_KR.html</td></tr>
      <tr><td>상세 리뷰 페이지</td><td>gz/results/review/manuscript_review_KR.html</td></tr>
      <tr><td>계획서 · 감사 기록</td><td>gz/plan/ · logs/AUDIT_LOG.md</td></tr>
      <tr><td>공개 저장소</td><td>github.com/vincentjshim-rgb/2C_ZGA</td></tr>
    </tbody>
  </table>
</div>

<footer>
  생성 2026-09-20 · 원고 <code>MANUSCRIPT_gz_v1_KR.md</code> / <code>MANUSCRIPT_gz_v1_EN.md</code> 기준 ·
  모든 수치는 <code>logs/AUDIT_LOG.md</code>로 추적된다 · 이 페이지는 설명용이며 투고본이 아니다.
</footer>

</div>
'''

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, 'w', encoding='utf-8') as fh:
    fh.write(HTML)
print('wrote', OUT, os.path.getsize(OUT) // 1024, 'KB')
