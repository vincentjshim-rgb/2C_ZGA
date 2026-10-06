#!/usr/bin/env python3
"""Build the general-reader explainer (public_KR.html), self-contained with embedded panels."""
import base64
import os

ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
B = f'{ANALYSIS}/gz'
ASSETS = f'{B}/figures/review_assets'
CROP = f'{ASSETS}/pubcrop'
OUT = f'{B}/results/review/public_KR.html'


def img(key):
    with open(f'{CROP}/{key}.jpg', 'rb') as fh:
        return 'data:image/jpeg;base64,' + base64.b64encode(fh.read()).decode()


def fig(key, cap, wide=False):
    return (f'<figure class="panel{" wide" if wide else ""}">'
            f'<img src="{img(key)}" alt="{cap}" loading="lazy">'
            f'<figcaption>{cap}</figcaption></figure>')


STEPS = [
    dict(
        n='1', tag='있는가', cls='s1',
        q='그 하락이 정말 있는가, 그리고 언제 일어나는가?',
        why='한 연구의 우연이라면 나머지 질문은 의미가 없다. 그래서 가장 먼저, 서로 다른 데이터에서 같은 일이 '
            '일어나는지부터 확인한다.',
        body=[
            ('html', fig('f1b_all', '그림 1b — 네 종의 시계 값. 점 하나가 배아 하나, 마름모가 단계 평균이다. '
                                    '연한 띠가 문헌에서 가져온 ZGA 구간.', wide=True)),
            ('p', '네 종(생쥐·소·돼지·토끼)의 난자부터 상실배까지 시계를 대봤다.'),
            ('p', '<b>생쥐</b>에서 가장 큰 하락은 초기→후기 2세포기에 있었다(−0.19). 앞에서 말한 '
                  '<b>분열이 없는 바로 그 구간</b>이다. 다른 생쥐 데이터 두 건에서도 같은 구간에서 내려갔다'
                  '(−0.14, −0.26). <b>세 번 재현됐다.</b>'),
            ('html', fig('f1c', '그림 1c — 두 번째 생쥐 데이터. 검은 점은 전체 유전자, 파란 점은 발생 중 크게 '
                                '변하는 유전자를 빼고 다시 잰 값. 유전자를 빼도 2세포기 하락이 남는다.')),
            ('p', '<b>소·돼지·토끼</b>에서는 가장 큰 하락이 문헌상 ZGA 구간 <b>밖</b>에 있었다. '
                  '사전에 정한 기준은 4종 중 1종에서만 충족됐다.'),
            ('honest', '이 음성 결과를 왜 첫 그림에 그대로 두는가',
             '미리 정한 규칙의 결과이기 때문이다. 맞으면 쓰고 틀리면 빼는 것은 조작이다. '
             '나중에 문헌을 뒤져보니 이유가 있었다 — <b>생쥐 구간만 분열이 없고 나머지 세 종의 구간에는 '
             '분열이 들어 있다.</b> 애초에 같은 검정이 아니었던 것이다. 이 사실도 논문에 적었다.'),
        ]),
    dict(
        n='2', tag='무엇의 합인가', cls='s2',
        q='그런데 그 "숫자 하나"는 대체 어떻게 생긴 걸까?',
        why='지금까지는 단계마다 숫자 하나만 봤다. 그 숫자가 유전자 1,839개의 덧셈이라면, 덧셈되는 '
            '<b>행렬 자체</b>를 보여주는 것이 맞다.',
        body=[
            ('html', fig('f1e', '그림 1e — 왼쪽: 가로가 발생 단계, 세로가 시계 유전자 1,839개. 가운데: 각 유전자가 '
                                '2세포기 하락에 얼마나 기여하는지. 오른쪽: 그것을 전부 가중합한 것 — 지금까지 봐온 '
                                '"시계 값"이다.', wide=True)),
            ('p', '유전자를 기여도 순으로 줄 세웠다. 한눈에 보이는 것은 <b>하락을 나르는 유전자가 맨 위와 맨 아래에 '
                  '몰려 있고, 가운데의 대다수는 거의 기여하지 않는다</b>는 점이다.'),
            ('p', '오른쪽 선이 왼쪽 행렬에서 계산된 것이다. 이 두 그림이 <b>같은 데이터</b>라는 사실이 다음 질문들의 '
                  '출발점이 된다.'),
        ]),
    dict(
        n='3', tag='청소 효과인가', cls='s3',
        q='그냥 엄마 RNA가 사라진 것 아닌가?',
        why='앞에서 말한 <b>함정 1(비율의 함정)</b>을 직접 검정한다. 이걸 통과하지 못하면 이 연구의 나머지는 '
            '의미가 없다.',
        body=[
            ('p', '방법은 이렇다. 접합자 쪽은 <b>아무것도 바꾸지 않고</b>, 난자에서 모계 유전자만 실제로 관찰된 '
                  '비율만큼 줄인 뒤 시계를 다시 매겼다. 즉 <b>"청소만 일어났다면 시계는 얼마나 움직였을까"</b>를 '
                  '계산한 것이다.'),
            ('html', fig('f1d', '그림 1d — 맨 오른쪽 "sim."이 청소만 모의한 결과다. 원래 하락(V0)과 비교해 보면 '
                                '상당 부분이 그것만으로 재현된다.')),
            ('num', ['원래 하락의 <b>0.70</b>이 청소만으로 재현 (생쥐)',
                     '<b>0.47</b> 재현 (소)']),
            ('p', '<b>절반쯤은 청소 효과가 맞다.</b> 논문은 이것을 먼저 인정한다.'),
            ('honest', '왜 스스로 약점을 먼저 말하는가',
             '어차피 심사자가 물어볼 것이다. 우리가 먼저 재면 <b>남은 절반</b>이 비로소 진짜 질문이 된다. '
             '숨기면 논문 전체가 이 한 가지로 무너지지만, 먼저 재면 그 위에 다음 질문을 쌓을 수 있다.'),
        ]),
    dict(
        n='4', tag='ZGA에 달렸나', cls='s4',
        q='그럼 남은 부분은 ZGA에 달려 있는가?',
        why='여기가 이 논문이 선행 연구와 갈라지는 지점이다. 지금까지의 배아 시계 연구는 모두 <b>정상 배아를 '
            '관찰</b>하기만 했다. ZGA를 실제로 막고 시계를 대본 사람이 없었다.',
        body=[
            ('html', fig('f2a', '그림 2a — 실험 설계. 색칠된 칸이 교란이 걸려 있는 단계, 점이 라이브러리를 딴 '
                                '단계다.', wide=True)),
            ('p', 'A485라는 약으로 ZGA를 막은 배아와 막지 않은 배아를 나란히 놓았다. 점 하나가 라이브러리 하나이고, '
                  '세 데이터셋의 세로축을 <b>같은 눈금</b>으로 맞췄기 때문에 하락의 깊이를 눈으로 비교할 수 있다.'),
            ('html', fig('f2b', '그림 2b — 빈 원이 초기 2세포기, 찬 원이 후기 2세포기. arm 위의 D가 그 arm의 하락, '
                                '아래 I가 대조군과의 차이다.', wide=True)),
            ('num', ['대조군 <b>−0.24</b> → A485 <b>−0.06</b>',
                     '차이 I = <b>+0.18</b> [0.11, 0.25]']),
            ('p', '하락이 4분의 1로 줄었다.'),
            ('honest', '그런데 사전 규칙은 충족되지 않았다',
             '모계 Brg1 제거도 같은 방향(+0.11)이지만 구간이 0을 포함한다. 가장 작은 데이터(Tardbp 제거)는 '
             '+0.04이고, 유전자 집합을 바꾸면 부호가 뒤집힌다. 규칙이 "판정 가능한 모든 데이터에서"를 요구했으므로 '
             '<b>전체로는 미충족</b>이다.<br><br>'
             '왜 그 데이터를 빼지 않는가: <b>결과를 보고 데이터를 빼는 것이 바로 사전 명시가 막으려는 행위다.</b> '
             '게다가 그 데이터가 아무것도 말하지 않는다는 사실 자체가 정보다 — 신뢰구간의 폭이 0.21인데 '
             '움직인 거리는 0.056이었다.'),
        ]),
    dict(
        n='5', tag='무엇으로 이루어졌나', cls='s5',
        q='그 숫자는 무엇으로 이루어져 있는가?',
        why='앞에서 말한 <b>함정 2(상쇄의 함정)</b>를 푼다. 그리고 여기가 이 논문의 심장이다.',
        body=[
            ('formula', '어떤 유전자의 기여 = (그 유전자의 시계 가중치) × (그 유전자의 발현 변화)'),
            ('p', '이 항들을 전부 더하면 보고된 시계 차이와 <b>소수점 아홉 자리까지</b> 같다. '
                  '추정이나 모델이 아니라 <b>항등식</b>이다 — 반박할 여지가 없다는 뜻이다.'),
            ('html', fig('f3a', '그림 3a — 유전자를 기여도 순으로 줄 세우고 왼쪽부터 차례로 더해 나간 그래프. '
                                '끝점이 논문이 보고하는 값이다.')),
            ('p', '<b>U자 모양이 맞다.</b> 음수부터 더하니 −1.40까지 내려갔다가, 양수 구간에서 올라와 '
                  '−0.24에서 끝난다. 그 끝점이 우리가 "시계 값"이라고 불러온 숫자다.'),
            ('num', ['내리는 힘 <b>−1.40</b>', '올리는 힘 <b>+1.17</b>', '남는 것 <b>−0.24</b> (17%)']),
            ('html', fig('f3b', '그림 3b — 하락에 가장 크게 기여한 20개 유전자(검정)와, 같은 유전자의 A485에서의 '
                                '기여(황토색). 오른쪽 숫자는 그 유전자의 발현 변화량.')),
            ('p', '<b>A485가 지운 것은 바로 그 지배적 유전자들이었다.</b> 상위 10개가 전체 차이의 63%를 설명한다.'),
            ('html', fig('f3c', '그림 3c — 연구실도, 교란 방법도, 라이브러리 종류도 다른 두 연구에서 유전자별 기여를 '
                                '맞대본 것.')),
            ('p', '두 독립 연구에서 <b>r = 0.74</b>로 일치한다. 이 분해가 우연이나 잡음이 아니라는 뜻이다.'),
        ]),
    dict(
        n='6', tag='되살리면?', cls='s6', key=True,
        q='ZGA를 되살리면 시계도 돌아오는가?',
        why='시계가 "재프로그래밍 상태"를 제대로 재는 도구라면, ZGA를 되살렸을 때 시계도 돌아와야 한다. '
            '<b>안 돌아온다면 시계가 그것을 읽지 못한다는 뜻이다.</b> 이 실험이 결정적인 이유다.',
        body=[
            ('p', 'A485로 막은 상태에서 <b>DUX</b>라는 인자를 넣어 접합자 전사를 되살렸다.'),
            ('html', fig('f3d', '그림 3d — 접합자 유전자와 2세포기 마커 유전자의 발현. 빨강이 높음, 파랑이 낮음. '
                                '오른쪽 구제 arm에서 붉은 블록이 돌아온 것이 보인다.', wide=True)),
            ('html', fig('f3e', '그림 3e — 접합자 유전자 점수를 숫자 하나로. 구제 arm에서 확실히 올라갔다.')),
            ('p', '접합자 유전자 발현이 <b>실제로 돌아왔다</b>(A485 대비 +1.67).'),
            ('html', fig('f3f', '그림 3f — 세 arm의 누적 기여 곡선을 <b>대조군 순서 그대로</b> 쌓은 것. '
                                '같은 유전자를 비교하기 위해서다. 점선이 대조군에서 값을 내리던 666개 유전자 지점.')),
            ('p', '점선 위치에서 내리는 기여를 읽으면 <b>−1.40 → −0.76 → −1.19</b>. 구제가 <b>85%</b>를 '
                  '되돌렸다.'),
            ('p', '그런데 곡선의 <b>끝점</b>을 보라. <b>−0.24 / −0.06 / −0.06</b>.'),
            ('punch', '구제했는데, 시계 값은 A485와 완전히 똑같다.'),
            ('p', '왜인가? 올리는 기여도 같이 커졌기 때문이다(대조군의 113%). 앞에서 든 회계 비유가 여기에 '
                  '그대로 들어맞는다.'),
            ('table', ['', '내리는 힘', '올리는 힘', '순 값'],
             [['A485', '−1.15', '+1.10', '<b>−0.06</b>'],
              ['A485 + DUX', '−1.38', '+1.32', '<b>−0.06</b>']],
             '순 값만 보면 두 arm을 구분할 수 없다. 그러나 안쪽의 규모는 20% 이상 다르다.'),
            ('p', '<b>이것이 이 논문이 말하려는 것이다.</b> 시계 값 하나로 "회복되었다 / 안 되었다"를 판정하면, '
                  '실제로 회복된 것을 놓친다.'),
        ]),
    dict(
        n='7', tag='다른 자로', cls='s7',
        q='시계를 쓰지 않는 자로 재도 같은 답이 나오는가?',
        why='여기까지의 모든 결론이 시계 하나에 의존한다. 완전히 다른 원리의 자가 하나 더 필요하다.',
        body=[
            ('aside', 'splicing(이어맞추기)이란',
             '유전자에서 막 만들어진 RNA에는 쓰지 않는 조각이 섞여 있다. 세포는 그것을 잘라내고 필요한 부분만 '
             '이어붙인다. 이어붙이는 방식은 한 가지가 아니어서 <b>같은 유전자에서 여러 버전의 RNA</b>가 나온다. '
             '배아 발생 중에 이 방식이 바뀐다.'),
            ('p', '여기에는 시계 가중치가 <b>한 번도</b> 등장하지 않는다. 완전히 독립된 측정이다.'),
            ('html', fig('f4c', '그림 4c — 각 라이브러리를 "대조군 초기 = 0, 대조군 후기 = 1" 축 위에 올린 '
                                '진행 점수.')),
            ('num', ['A485 <b>0.56</b>', 'Tardbp 제거 <b>0.58</b>', 'Brg1 제거 <b>0.55</b>']),
            ('p', '세 교란 모두 진행이 절반 수준으로 줄었고, 구제하면 대조군 쪽으로 돌아온다. '
                  '<b>그리고 이 지표는 사전 규칙을 충족했다 — 시계 지표는 못 했다.</b>'),
            ('honest', 'splicing이 더 좋은 지표라는 뜻이 아니다',
             '두 규칙의 난이도가 달랐다(시계 쪽이 더 빡셌다). 그리고 splicing 점수는 수백 개 이벤트를 '
             '<b>같은 방향으로 맞춰</b> 평균 내므로 상쇄가 일어나지 않는다. 시계는 부호가 섞인 합이라 상쇄된다.<br><br>'
             '즉 이 대비 자체가 <b>"덧셈식은 상쇄된다"</b>는 이 논문의 메시지를 뒷받침한다. '
             'splicing 쪽에도 약점이 있었다 — 표본 내 채점이 크기를 부풀려서, 교차 검증 후 방향만 쓰기로 했다.'),
            ('html', fig('f4e', '그림 4e — 어떤 종류의 이어맞추기가 바뀌는지. 왼쪽은 활성화된 비율, 오른쪽은 '
                                '그중 증가한 비율.', wide=True)),
            ('p', '흥미로운 점: 활성화된 이벤트의 <b>43~47%만 증가</b>했다. 즉 "splicing 활성화 = 무조건 증가"가 '
                  '아니라 양방향이다.'),
        ]),
    dict(
        n='8', tag='read로 확인', cls='s8',
        q='추정이 아니라 실제 read로 세어도 같은가?',
        why='splicing을 <b>추정치</b>로 쟀다는 것이 이 논문의 약점 중 하나였다. 그런데 이건 '
            '<b>직접 검증할 수 있는</b> 약점이다. 그래서 검증했다.',
        body=[
            ('p', '23개 라이브러리를 유전체에 다시 정렬해서, 이어붙인 자리를 지나가는 <b>실제 read 개수</b>를 셌다. '
                  '다른 프로그램, 다른 통계, 다른 정량법이다.'),
            ('html', fig('s7a', '보충 그림 S7a — 두 방법으로 잰 PSI를 맞대본 밀도 지도. 대각선 능선이 일치를 '
                                '보여준다.')),
            ('num', ['58,958쌍에서 <b>ρ = 0.78</b>', '교란·구제의 방향도 재현']),
            ('html', fig('s7b', '보충 그림 S7b — 실제 read를 직접 본 것. 산이 read가 쌓인 곳, 호 위의 숫자가 '
                                '그 자리를 건너간 read 개수다.')),
            ('table', ['', '이어붙임 1', '이어붙임 2'],
             [['대조군', '414', '514'], ['A485', '<b>74</b>', '<b>124</b>'], ['A485 + DUX', '318', '415']],
             '평균이나 통계가 아니라 read 개수로 보인다. A485에서 무너지고 구제에서 돌아온다.'),
        ]),
    dict(
        n='9', tag='사람도 그런가', cls='s9', key=True,
        q='사람 배아에서도 같은 일이 일어나는가?',
        why='여기까지는 전부 <b>생쥐 한 종, 2세포기 한 구간</b>이다. "거의 상쇄되고 남은 나머지"가 이 구간만의 '
            '특징인지, 아니면 시계라는 도구가 배아를 읽을 때 으레 그렇게 되는 것인지 아직 모른다. '
            '구분하려면 다른 종, 다른 구간에서 같은 계산을 해봐야 한다.',
        body=[
            ('p', '사람 배아 공개 데이터를 가져와 배아 20개로 묶고, <b>같은 시계를 같은 방식으로</b> 댔다. '
                  '볼 구간은 결과를 보기 전에 <b>8세포기 → 상실배</b>로 정해 잠갔다. 선행연구가 사람에서 하락을 '
                  '보고한 구간이 거기이기 때문이다.'),
            ('html', fig('f5a', '그림 5a — 사람 배아의 시계 값. 점 하나가 배아 하나, 마름모가 단계 평균이다. '
                                '연한 띠가 미리 정해 둔 구간(8세포기→상실배).')),
            ('p', '그 구간에서 시계가 <b>−0.10</b> 내려갔다. 사람 궤적에서 가장 큰 하락이다.'),
            ('html', fig('f5b', '그림 5b — 사람의 같은 구간에서, 유전자 1,839개의 기여를 순서대로 더해 나간 것. '
                                '생쥐의 그림 3b와 같은 계산이다.')),
            ('num', ['내리는 유전자 696개가 <b>−1.24</b>', '올리는 유전자 674개가 <b>+1.15</b>',
                     '보고된 값 −0.10은 그 <b>8%</b>']),
            ('punch', '생쥐 17%, 사람 8%. 남는 비율까지 비슷하다 — <b>구조가 같다.</b>'),
            ('p', '더 중요한 것은 조건이 같지 않다는 점이다. 사람의 이 구간에는 <b>세포분열이 들어 있다.</b> '
                  '생쥐 무대의 특징(분열이 없다)이 없는데도 같은 모양이 나왔다. 그렇다면 이 "거의 상쇄"는 '
                  '2세포기라는 특정 시기의 성질이 아니라, <b>전사체가 통째로 갈릴 때 가중치 덧셈식이 보이는 '
                  '성질</b>에 가깝다.'),
            ('html', fig('f5c', '그림 5c — 유전자별 기여도를 사람(세로)과 생쥐(가로)로 맞대본 것. '
                                '점이 대각선에 몰리면 같은 유전자가 같은 일을 한다는 뜻인데, 그렇지 않다.')),
            ('honest', '그런데 유전자 명단은 전혀 달랐다',
             '두 종에서 각 유전자가 낸 기여도는 거의 무관했다(<b>r = 0.07</b>). 같은 구간을 본 생쥐 연구 '
             '<b>두 건끼리는 r = 0.74</b>였으니, 계산이 흔들려서 생긴 차이가 아니다. '
             '정리하면 — <b>“거의 상쇄된다”는 성질은 종과 구간을 넘어 반복되지만, 그것을 만드는 유전자 명단은 '
             '반복되지 않는다.</b> 그래서 기여도 상위 유전자를 "배아 노화의 주역"처럼 읽으면 안 된다. '
             '그건 이 시계가 이 구간을 읽는 방식이지, 유전자의 순위표가 아니다.'),
        ]),
]

LIMITS = [
    ('사람 쪽 비교는 <b>배아 다섯 개</b>다', '8세포기 3개, 상실배 2개. 하락 −0.10의 구간 '
     '추정은 0을 포함한다. 또 사람 유전자를 생쥐 유전자에 대응시켜야 했는데, 그 대응표는 시계를 만든 쪽이 '
     '함께 배포한 것을 그대로 썼다.'),
    ('교란 그룹당 배아 라이브러리가 <b>2~4개</b>다', '공개 데이터에 그 이상이 없다. 그래서 작은 p 값이 아니라 '
     '효과 크기·구간·데이터 간 방향 일치로 논증한다.'),
    ('약물과 유전자 제거는 <b>ZGA만</b> 건드리지 않는다', 'A485는 p300/CBP라는 효소를 막는데, 그 효소는 ZGA 말고도 '
     '많은 일을 한다. 그래서 "ZGA가 원인이다"라고 쓰지 않고 "ZGA를 막으면 약해진다"까지만 쓴다.'),
    ('핵 이식(복제) 배아 결과는 <b>설명하지 못한다</b>', '복제 배아가 체외수정 배아보다 낮게 나왔는데, 이는 사전 '
     '예상과 반대다. 이유를 모르므로 그대로 적고 부호를 열어 두었다.'),
    ('분해·구제·교차검증은 <b>사후 분석</b>이다', '각각 실행 전에 계획을 썼지만 판정을 본 뒤였다. 해당 절마다 '
     '"(사후)"로 표시했다.'),
    ('시계는 <b>성체 조직</b>으로 훈련됐다', '배아에서의 판독은 그 훈련을 반영한다. 소·돼지·토끼 값은 훈련 종 '
     '밖이라 상동유전자 대응에 추가로 기댄다.'),
]

GLOSSARY = [
    ('노화 시계', '분자 프로파일을 읽어 생물학적 나이를 추정하는 공식. 유전자마다 가중치가 붙은 덧셈식이다.'),
    ('ZGA (접합자 유전체 활성화)', '배아가 자기 유전체를 처음 켜는 사건. 생쥐에서는 약한 물결과 강한 물결 두 번에 '
     '나눠 일어난다.'),
    ('2세포기', '수정란이 한 번 분열해 세포가 두 개가 된 시기. 초기와 후기 사이에 <b>분열이 없다</b>.'),
    ('모계 RNA', '엄마가 난자에 미리 넣어둔 RNA. 배아 초기에 대량 분해된다.'),
    ('DUX', '초기 배아 유전자를 켜는 전사 인자. 이 연구에서는 막힌 ZGA를 되살리는 데 썼다.'),
    ('A485', 'p300/CBP라는 효소의 활성을 막는 약. 이 연구에서 ZGA 차단에 사용했다.'),
    ('splicing (이어맞추기)', 'RNA에서 불필요한 부분을 잘라내고 필요한 부분을 이어붙이는 과정.'),
    ('PSI', '어떤 조각이 최종 RNA에 포함된 비율. splicing 방식을 숫자로 나타낸 것.'),
    ('상실배 (morula)', '8세포기 다음 단계. 세포가 오디처럼 뭉친 모양이라 붙은 이름이다. 사람에서 시계가 내려간 구간의 끝이다.'),
]

WHOLE = [('whole_F1', '그림 1'), ('whole_F2', '그림 2'), ('whole_F3', '그림 3'),
         ('whole_F4', '그림 4'), ('whole_F5', '그림 5'), ('whole_S7', '보충 그림 S7')]

# ------------------------------------------------------------------------ render


def block(b):
    k = b[0]
    if k == 'p':
        return f'<p>{b[1]}</p>'
    if k == 'html':
        return b[1]
    if k == 'formula':
        return f'<div class="formula">{b[1]}</div>'
    if k == 'punch':
        return f'<p class="punch">{b[1]}</p>'
    if k == 'num':
        return '<ul class="nums">' + ''.join(f'<li>{x}</li>' for x in b[1]) + '</ul>'
    if k == 'honest':
        return f'<div class="honest"><h4>{b[1]}</h4><p>{b[2]}</p></div>'
    if k == 'aside':
        return f'<div class="aside"><h4>{b[1]}</h4><p>{b[2]}</p></div>'
    if k == 'table':
        head = ''.join(f'<th>{h}</th>' for h in b[1])
        rows = ''.join('<tr>' + ''.join(f'<td>{c}</td>' for c in r) + '</tr>' for r in b[2])
        cap = f'<figcaption>{b[3]}</figcaption>' if len(b) > 3 else ''
        return f'<div class="tbl"><table><thead><tr>{head}</tr></thead><tbody>{rows}</tbody></table>{cap}</div>'
    raise ValueError(k)


steps_html = '\n'.join(
    f'''<section class="step {s['cls']}{' key' if s.get('key') else ''}" id="step{s['n']}">
  <div class="stephead">
    <span class="stepnum">{s['n']}</span>
    <span class="steptag">{s['tag']}</span>
  </div>
  <h2>{s['q']}</h2>
  <div class="why"><span class="whylab">왜 이 질문인가</span><p>{s['why']}</p></div>
  {''.join(block(b) for b in s['body'])}
</section>''' for s in STEPS)

chain_html = '\n'.join(
    f'<li><span class="cn">{s["n"]}</span><span class="cq">{s["q"]}</span></li>' for s in STEPS)

limits_html = '\n'.join(f'<div class="lim"><h4>{a}</h4><p>{b}</p></div>' for a, b in LIMITS)
gloss_html = '\n'.join(f'<div class="gl"><dt>{a}</dt><dd>{b}</dd></div>' for a, b in GLOSSARY)
whole_html = '\n'.join(
    f'<figure class="whole"><img src="{img(k)}" alt="{c}" loading="lazy"><figcaption>{c}</figcaption></figure>'
    for k, c in WHOLE)

HTML = f'''<title>배아는 정말 젊어지는가</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Nanum+Myeongjo:wght@400;700;800&family=Noto+Sans+KR:wght@300;400;500;700&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
:root {{
  --paper:#FAF9F6; --card:#FFFFFF; --ink:#1C1B19; --ink2:#3D3A35; --muted:#6F6B64;
  --rule:#E6E1D8; --rule2:#D3CCBF;
  --amber:#A9700C; --teal:#12786E; --down:#A33B2B; --up:#2E5B87; --violet:#5B4A86;
  --band:#F2E9D4;
  --shadow:0 1px 2px rgba(28,27,25,.04), 0 10px 30px -20px rgba(28,27,25,.35);
}}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{
    --paper:#141310; --card:#1C1A16; --ink:#EFEBE2; --ink2:#D3CDC2; --muted:#9B958A;
    --rule:#2D2A24; --rule2:#3C382F;
    --amber:#D89A3E; --teal:#4FB8AC; --down:#D57B6B; --up:#7FA9D0; --violet:#A797D4;
    --band:#2E2919;
    --shadow:0 1px 2px rgba(0,0,0,.4), 0 10px 30px -20px rgba(0,0,0,.9);
  }}
}}
:root[data-theme="dark"] {{
  --paper:#141310; --card:#1C1A16; --ink:#EFEBE2; --ink2:#D3CDC2; --muted:#9B958A;
  --rule:#2D2A24; --rule2:#3C382F;
  --amber:#D89A3E; --teal:#4FB8AC; --down:#D57B6B; --up:#7FA9D0; --violet:#A797D4;
  --band:#2E2919;
  --shadow:0 1px 2px rgba(0,0,0,.4), 0 10px 30px -20px rgba(0,0,0,.9);
}}
* {{ box-sizing:border-box; }}
body {{
  margin:0; background:var(--paper); color:var(--ink);
  font-family:'Noto Sans KR', -apple-system, 'Apple SD Gothic Neo', 'Malgun Gothic', sans-serif;
  font-size:16.5px; line-height:1.95; word-break:keep-all; -webkit-font-smoothing:antialiased;
}}
.wrap {{ max-width:1120px; margin:0 auto; padding-left:20px; padding-right:20px; padding-block:0 110px; }}
h1,h2,h3,h4 {{ font-family:'Nanum Myeongjo','Noto Serif KR',Georgia,serif; text-wrap:balance; margin:0; }}
p {{ margin:0 0 1.05em; }}
b {{ font-weight:700; color:var(--ink); }}
.measure {{ max-width:760px; }}

header.top {{ position:sticky; top:0; z-index:30; background:color-mix(in srgb, var(--paper) 93%, transparent);
  backdrop-filter:blur(10px); border-bottom:1px solid var(--rule); }}
.topin {{ max-width:1120px; margin:0 auto; padding:11px 20px; display:flex; gap:16px; align-items:baseline;
  flex-wrap:wrap; }}
.brand {{ font-family:'Nanum Myeongjo',serif; font-weight:800; font-size:14.5px; white-space:nowrap; }}
nav {{ display:flex; gap:3px; flex-wrap:wrap; }}
nav a {{ font-size:12.5px; color:var(--muted); text-decoration:none; padding:3px 9px; border-radius:999px;
  border:1px solid transparent; }}
nav a:hover {{ color:var(--ink); border-color:var(--rule2); }}

.hero {{ padding-block:68px 46px; border-bottom:1px solid var(--rule); }}
.kicker {{ font-family:'IBM Plex Mono',monospace; font-size:11.5px; letter-spacing:.14em; text-transform:uppercase;
  color:var(--amber); margin-bottom:20px; }}
.hero h1 {{ font-size:clamp(32px,5.4vw,50px); line-height:1.3; margin-bottom:24px; }}
.hero .lede {{ font-size:18px; color:var(--ink2); max-width:720px; }}
.oneline {{ margin-top:32px; padding:22px 26px; border-left:3px solid var(--amber); background:var(--card);
  border-radius:0 12px 12px 0; box-shadow:var(--shadow); }}
.oneline .l {{ font-family:'IBM Plex Mono',monospace; font-size:11px; letter-spacing:.11em; color:var(--muted);
  text-transform:uppercase; display:block; margin-bottom:8px; }}
.oneline p {{ margin:0; font-size:17.5px; }}

h2.sec {{ font-size:27px; margin:76px 0 8px; }}
h2.sec + .sub {{ color:var(--muted); font-size:14.5px; margin-bottom:30px; max-width:760px; }}

.bg {{ background:var(--card); border:1px solid var(--rule); border-radius:14px; padding:28px 30px; margin:22px 0;
  box-shadow:var(--shadow); }}
.bg h3 {{ font-size:20px; margin-bottom:14px; }}
.bg .eyebrow {{ font-family:'IBM Plex Mono',monospace; font-size:11px; letter-spacing:.12em;
  text-transform:uppercase; color:var(--amber); display:block; margin-bottom:8px; }}
.bg p {{ font-size:15.8px; color:var(--ink2); }}
.bg ol {{ margin:0 0 1em; padding-left:22px; }}
.bg ol li {{ font-size:15.8px; color:var(--ink2); margin-bottom:5px; }}

.formula {{ font-family:'IBM Plex Mono',monospace; font-size:14px; background:var(--band); color:var(--ink);
  border-radius:10px; padding:16px 20px; margin:18px 0; line-height:1.7; overflow-x:auto; }}

.analogy {{ border:1px dashed var(--rule2); border-radius:12px; padding:20px 24px; margin:18px 0;
  background:color-mix(in srgb, var(--band) 45%, transparent); }}
.analogy h4 {{ font-size:14px; color:var(--amber); margin-bottom:9px; font-family:'IBM Plex Mono',monospace;
  letter-spacing:.06em; }}
.analogy p {{ margin:0 0 .6em; font-size:15.5px; color:var(--ink2); }}
.analogy p:last-child {{ margin-bottom:0; }}

.trap {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(310px,1fr)); gap:20px; }}
.trap .t {{ background:var(--card); border:1px solid var(--rule); border-top:3px solid var(--down);
  border-radius:12px; padding:24px 24px 20px; }}
.trap .t:last-child {{ border-top-color:var(--up); }}
.trap .t .n {{ font-family:'IBM Plex Mono',monospace; font-size:11px; letter-spacing:.1em; color:var(--muted); }}
.trap .t h4 {{ font-size:18px; margin:7px 0 11px; }}
.trap .t p {{ font-size:15.2px; color:var(--ink2); margin-bottom:.7em; }}

.chain {{ list-style:none; margin:0; padding:0; border-top:1px solid var(--rule); }}
.chain li {{ display:grid; grid-template-columns:40px 1fr; gap:14px; padding:14px 0;
  border-bottom:1px solid var(--rule); align-items:baseline; }}
.cn {{ font-family:'IBM Plex Mono',monospace; font-size:13px; color:var(--amber); }}
.cq {{ font-size:16px; color:var(--ink2); }}

.step {{ margin:60px 0; padding-top:34px; border-top:1px solid var(--rule); }}
.step.key {{ border-top:2px solid var(--amber); }}
.stephead {{ display:flex; align-items:center; gap:12px; margin-bottom:12px; }}
.stepnum {{ font-family:'IBM Plex Mono',monospace; font-size:13px; width:30px; height:30px; border-radius:50%;
  display:grid; place-items:center; background:var(--card); border:1px solid var(--rule2); color:var(--amber); }}
.steptag {{ font-size:12px; padding:3px 11px; border-radius:999px; border:1px solid var(--rule2);
  color:var(--ink2); background:var(--card); }}
.step h2 {{ font-size:clamp(23px,3vw,30px); line-height:1.45; margin-bottom:20px; max-width:860px; }}
.step p, .step ul {{ max-width:760px; }}

.why {{ max-width:760px; margin:0 0 26px; padding:16px 0 16px 20px; border-left:3px solid var(--amber); }}
.whylab {{ font-family:'IBM Plex Mono',monospace; font-size:11px; letter-spacing:.1em; text-transform:uppercase;
  color:var(--amber); display:block; margin-bottom:7px; }}
.why p {{ margin:0; font-size:16px; color:var(--ink2); }}

.panel {{ margin:26px 0; max-width:760px; }}
.panel.wide {{ max-width:1080px; }}
.panel img {{ display:block; width:100%; height:auto; background:#fff; border:1px solid var(--rule);
  border-radius:8px; padding:10px; }}
.panel figcaption {{ font-size:13.5px; color:var(--muted); margin-top:10px; line-height:1.7; }}

.punch {{ font-family:'Nanum Myeongjo',serif; font-size:clamp(21px,2.8vw,27px); line-height:1.55;
  color:var(--ink); border-top:1px solid var(--rule2); border-bottom:1px solid var(--rule2);
  padding:24px 0; margin:30px 0; max-width:760px; }}

ul.nums {{ list-style:none; margin:16px 0 20px; padding:0; display:flex; flex-wrap:wrap; gap:8px; }}
ul.nums li {{ font-family:'IBM Plex Mono',monospace; font-size:13px; color:var(--ink2); background:var(--card);
  border:1px solid var(--rule); border-radius:7px; padding:7px 13px; }}
ul.nums li b {{ font-family:'IBM Plex Mono',monospace; }}

.honest {{ max-width:760px; margin:24px 0; padding:20px 24px; border-radius:0 12px 12px 0;
  border-left:3px solid var(--violet); background:color-mix(in srgb, var(--violet) 7%, transparent); }}
.honest h4 {{ font-size:16px; margin-bottom:9px; }}
.honest p {{ margin:0; font-size:15.2px; color:var(--ink2); }}

.aside {{ max-width:760px; margin:24px 0; padding:20px 24px; border-radius:12px; background:var(--card);
  border:1px solid var(--rule); }}
.aside h4 {{ font-size:16px; margin-bottom:9px; color:var(--teal); }}
.aside p {{ margin:0; font-size:15.2px; color:var(--ink2); }}

.tbl {{ max-width:760px; margin:22px 0; }}
.tbl table {{ border-collapse:collapse; width:100%; font-size:14.5px; }}
.tbl th, .tbl td {{ text-align:right; padding:11px 14px; border-bottom:1px solid var(--rule); }}
.tbl th:first-child, .tbl td:first-child {{ text-align:left; }}
.tbl th {{ font-family:'IBM Plex Mono',monospace; font-size:11px; letter-spacing:.07em; text-transform:uppercase;
  color:var(--muted); font-weight:500; }}
.tbl td {{ font-family:'IBM Plex Mono',monospace; color:var(--ink2); }}
.tbl td:first-child {{ font-family:'Noto Sans KR',sans-serif; }}
.tbl figcaption {{ font-size:13.5px; color:var(--muted); margin-top:10px; }}

.limgrid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(300px,1fr)); gap:18px; }}
.lim {{ background:var(--card); border:1px solid var(--rule); border-left:3px solid var(--down);
  border-radius:0 12px 12px 0; padding:20px 22px; }}
.lim h4 {{ font-size:15.5px; margin-bottom:8px; }}
.lim p {{ margin:0; font-size:14.5px; color:var(--ink2); }}

.glgrid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(320px,1fr)); gap:0 28px; }}
.gl {{ padding:14px 0; border-bottom:1px solid var(--rule); }}
.gl dt {{ font-weight:700; font-size:15.5px; margin-bottom:4px; }}
.gl dd {{ margin:0; font-size:14.5px; color:var(--ink2); }}

.wholegrid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(300px,1fr)); gap:24px; }}
.whole img {{ display:block; width:100%; height:auto; background:#fff; border:1px solid var(--rule);
  border-radius:8px; padding:8px; }}
.whole figcaption {{ font-size:13px; color:var(--muted); margin-top:8px; }}
.whole {{ margin:0; }}

footer {{ margin-top:80px; padding-top:24px; border-top:1px solid var(--rule); color:var(--muted);
  font-size:13px; max-width:760px; }}
footer a {{ color:var(--ink2); }}

@media (max-width:700px) {{
  body {{ font-size:16px; }}
  .step {{ margin:44px 0; }}
  .bg {{ padding:22px 18px; }}
  .tbl th, .tbl td {{ padding:9px 8px; }}
}}
</style>

<header class="top">
  <div class="topin">
    <span class="brand">배아는 정말 젊어지는가</span>
    <nav>
      <a href="#bg">배경</a>
      <a href="#trap">왜 어려운가</a>
      <a href="#chain">전체 흐름</a>
      <a href="#step1">1</a><a href="#step2">2</a><a href="#step3">3</a><a href="#step4">4</a>
      <a href="#step5">5</a><a href="#step6">6</a><a href="#step7">7</a><a href="#step8">8</a>
      <a href="#end">결론</a>
      <a href="#limits">한계</a>
      <a href="#gloss">용어</a>
    </nav>
  </div>
</header>

<div class="wrap">

<div class="hero measure">
  <div class="kicker">생쥐·사람 배아 · 분자 나이 · 2026</div>
  <h1>배아는 정말<br>“젊어지는가”</h1>
  <p class="lede">
    분자 나이를 재는 숫자 하나를 배아에 대보면 값이 <b>내려간다</b>. 나이가 거꾸로 간다는 뜻이다.
    이 연구는 그 숫자를 유전자 단위로 뜯어보고, 배아의 유전체를 실제로 막아보고, 사람 배아에서 다시 확인했다.
  </p>
  <div class="oneline">
    <span class="l">결론을 먼저</span>
    <p>하락은 진짜다. 그러나 그 숫자는 <b>서로 반대로 미는 두 힘이 거의 상쇄되고 남은 나머지</b>였다.
       생쥐에서도 사람에서도 그랬고, 그래서 안쪽에서 많은 것이 변해도 <b>숫자는 가만히 있을 수 있었다</b>.</p>
  </div>
</div>

<h2 class="sec" id="bg">먼저 알아야 할 두 가지</h2>
<p class="sub">이 연구를 이해하는 데 필요한 배경은 딱 두 가지다 — 노화 시계가 무엇인지, 그리고 배아의 첫 이틀에
   무슨 일이 일어나는지.</p>

<div class="bg measure">
  <span class="eyebrow">배경 1</span>
  <h3>노화 시계란 무엇인가</h3>
  <p>주민등록증의 나이와 몸의 나이는 다를 수 있다. 같은 70세라도 어떤 사람의 혈관은 60세처럼, 어떤 사람은
     80세처럼 보인다. 과학자들은 이 “몸의 나이”를 숫자로 재려고 <b>노화 시계</b>를 만들었다.</p>
  <p>만드는 방법은 의외로 단순하다.</p>
  <ol>
    <li>나이를 아는 시료를 수천 개 모은다</li>
    <li>각 시료에서 유전자 수만 개의 활동량을 잰다</li>
    <li>나이를 가장 잘 맞히는 <b>가중치 조합</b>을 컴퓨터로 찾는다</li>
  </ol>
  <p>결과물은 공식 하나다.</p>
  <div class="formula">예상 나이 = (유전자 A의 양 × 가중치 a) + (유전자 B의 양 × 가중치 b) + …</div>
  <p>이 연구가 쓴 시계는 생쥐·쥐·원숭이·사람의 <b>성체 조직 1만 1천 건 이상</b>으로 훈련된 것으로,
     유전자 <b>1,839개</b>에 가중치가 붙어 있다. 중요한 점 — <b>우리가 만든 게 아니라 이미 발표된 시계를
     그대로 가져다 썼다.</b> 이 논문의 목적은 새 시계를 만드는 것이 아니라, 기존 시계가 배아에서 무엇을
     재고 있는지 뜯어보는 것이기 때문이다.</p>
  <div class="analogy">
    <h4>왜 사람들이 놀랐나</h4>
    <p>2021년 이후 여러 연구가 이 시계를 초기 배아에 대보니 값이 <b>내려갔다</b>. “생애의 원점(ground zero)”
       이라는 말까지 나왔다.</p>
    <p>그럴싸한 이유가 있다. 부모는 <b>늙은 세포</b>에서 아이를 만든다. 그런데 아이는 <b>0살</b>부터 시작한다.
       그렇다면 어딘가에서 나이가 초기화되어야 한다. 배아가 그 지점이라면 말이 된다.</p>
  </div>
</div>

<div class="bg measure">
  <span class="eyebrow">배경 2</span>
  <h3>배아의 첫 이틀</h3>
  <p>난자는 비어 있는 세포가 아니다. 엄마가 미리 만들어 넣어준 RNA와 단백질이 <b>가득</b> 들어 있다.
     도시락을 싸서 보낸 셈이다.</p>
  <p>수정 직후 배아는 자기 유전체를 거의 읽지 않는다. 엄마 도시락으로 버틴다. 그러다 배아가
     <b>자기 유전체를 켠다</b>. 이것을 <b>접합자 유전체 활성화(ZGA)</b>라고 한다. 발생의 첫 번째 큰 사건이다.</p>
  <p>생쥐에서는 두 번에 나눠 일어난다. <b>약한 물결</b>(수정란~초기 2세포기)에 이어
     <b>강한 물결</b>(후기 2세포기)이 온다. 동시에 엄마 도시락은 대량으로 분해된다 —
     <b>들어오는 것과 나가는 것이 한꺼번에</b> 일어나는 시기다.</p>
  <div class="analogy">
    <h4>이 연구의 무대</h4>
    <p><b>초기 2세포기와 후기 2세포기 사이에는 세포분열이 없다.</b> 세포는 두 개 그대로이고 약 8시간이
       흐를 뿐이다.</p>
    <p>왜 중요한가 — 세포를 측정할 때 가장 흔한 반론이 “세포 구성이 바뀌어서 숫자가 바뀐 것 아니냐”이다.
       분열이 없으면 <b>같은 두 세포를 시간만 두고 본 것</b>이므로 그 반론이 성립하지 않는다.</p>
    <p>그리고 이건 생쥐의 행운이다. 소·돼지·토끼의 해당 구간에는 분열이 들어 있다.</p>
  </div>
</div>

<h2 class="sec" id="trap">왜 이 질문이 어려운가</h2>
<p class="sub">시계가 “가중치가 붙은 덧셈식”이라는 사실에서 두 가지 문제가 나온다.
   둘 다 정상 배아를 <b>관찰만 해서는 절대 풀 수 없다</b>.</p>

<div class="trap">
  <div class="t">
    <div class="n">함정 1</div>
    <h4>비율의 함정</h4>
    <p>유전자 활동량은 보통 “전체 RNA 중 몇 %”로 잰다. 그런데 배아 초기에는 엄마 RNA가 대량으로 사라진다.
       그러면 <b>아무 일도 하지 않은 유전자의 비율이 저절로 올라간다.</b></p>
    <div class="analogy">
      <p>교실에 100명이 있고 그중 남학생이 20명이다. 여학생 50명이 나가면 남학생 비율은 20%에서 40%가 된다.
         <b>남학생은 한 명도 늘지 않았는데.</b></p>
    </div>
    <p>시계 값이 움직였다고 해서 세포 상태가 변한 것은 아닐 수 있다.</p>
  </div>
  <div class="t">
    <div class="n">함정 2</div>
    <h4>상쇄의 함정</h4>
    <p>덧셈식에는 <b>음수와 양수가 섞여 있다.</b> 어떤 유전자는 값을 내리고 어떤 유전자는 올린다.
       둘 다 아주 크게 움직여도 <b>합은 거의 움직이지 않을 수 있다.</b></p>
    <div class="analogy">
      <p>A사: 매출 11억, 비용 11.5억 → 순이익 −0.5억<br>
         B사: 매출 13.2억, 비용 13.8억 → 순이익 −0.6억</p>
      <p>순이익만 보면 두 회사가 비슷하다. 그러나 <b>사업 규모는 20% 넘게 다르다.</b></p>
    </div>
    <p>시계 값은 “순이익”이다. 이 연구는 <b>매출과 비용을 따로 보기로 했다.</b></p>
  </div>
</div>

<div class="bg measure" style="margin-top:34px">
  <span class="eyebrow">이 연구의 방식</span>
  <h3>규칙을 먼저 쓰고 잠근다</h3>
  <p>결과를 본 뒤에 기준을 정하면 무엇이든 “성공”으로 만들 수 있다. 그래서 이 연구는 분석마다
     <b>질문·데이터·지표·판정 기준을 파일로 먼저 써서 잠갔다.</b> 잠근 날짜와 파일 지문이 공개 기록에 남아 있다.</p>
  <p>그 결과 <b>기준을 충족하지 못한 항목도 그대로 보고한다.</b> 아래에서 “사전 기준 미충족”이라는 말이 몇 번
     나오는데, 그건 실패를 숨기지 않았다는 뜻이다.</p>
</div>

<h2 class="sec" id="chain">전체 흐름</h2>
<p class="sub">아홉 개의 질문이 차례로 이어진다. 각 질문은 앞 질문의 답이 남긴 의심에서 나온다.</p>
<ul class="chain">{chain_html}</ul>

{steps_html}

<h2 class="sec" id="end">그래서 무엇을 알게 됐나</h2>
<div class="measure">
  <p>“배아에서 분자 나이가 내려간다”는 측정은 <b>진짜다</b>. 이 연구도 세 데이터에서 재현했다.</p>
  <p>그러나 그 숫자는—</p>
  <ol style="padding-left:22px">
    <li><b>절반쯤은</b> 엄마 RNA가 치워지는 효과다</li>
    <li>나머지는 ZGA를 막으면 <b>약해진다</b></li>
    <li>그리고 숫자 자체가 <b>−1.40과 +1.17이 상쇄되고 남은 17%의 나머지</b>다</li>
    <li>사람 배아의 다른 구간에서도 <b>같은 모양</b>(−1.24와 +1.15의 8%)이 나왔다 —
        단, <b>같은 유전자는 아니었다</b>(r = 0.07)</li>
  </ol>
  <p>그래서 안쪽에서 많은 것이 변해도 <b>숫자는 가만히 있을 수 있다</b>. DUX 구제 실험이 그 증거다.</p>
  <p>마지막 항목이 뜻하는 바는 이렇다. <b>“거의 상쇄된다”는 것은 배아의 성질이라기보다 시계라는 도구의
     성질에 가깝고</b>, 반대로 <b>기여도가 큰 유전자 명단은 그 시계가 그 구간을 읽는 방식일 뿐</b>이어서
     종을 옮기면 남지 않는다.</p>
  <div class="honest">
    <h4>실용적인 의미</h4>
    <p>회춘이나 노화 역전을 주장하는 연구가 시계 값 하나를 근거로 들 때, 그 값이 <b>무엇으로 이루어졌는지</b>
       함께 보아야 한다. 이 논문은 그렇게 보는 방법(유전자별 분해)과, 그렇게 봐야 하는 이유(구제 실험)를
       같이 보여준다.</p>
  </div>
</div>

<h2 class="sec" id="limits">정직하게 남긴 것</h2>
<p class="sub">논문이 스스로 적어 둔 한계들이다. 숨기면 심사에서 무너지고, 먼저 적으면 결론의 범위가 분명해진다.</p>
<div class="limgrid">{limits_html}</div>

<h2 class="sec" id="gloss">용어</h2>
<p class="sub">본문에 나온 순서대로.</p>
<div class="glgrid">{gloss_html}</div>

<h2 class="sec" id="whole">원본 그림</h2>
<p class="sub">위에서 조각으로 본 그림들의 전체 모습. 논문에 실리는 형태 그대로다.</p>
<div class="wholegrid">{whole_html}</div>

<footer>
  배아 전사체 나이 하락에 관한 연구의 일반 독자용 설명 · 2026-09-29 갱신<br>
  모든 그림은 논문에 실리는 것과 같은 그림이며, 수치는 공개 데이터에서 계산된 것이다.
  분석 코드와 사전 계획서, 감사 기록은 <a href="https://github.com/vincentjshim-rgb/2C_ZGA">github.com/vincentjshim-rgb/2C_ZGA</a>에 공개되어 있다.
</footer>

</div>
'''

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, 'w', encoding='utf-8') as fh:
    fh.write(HTML)
print('wrote', OUT, os.path.getsize(OUT) // 1024, 'KB')
