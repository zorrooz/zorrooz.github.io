<template>
  <div class="container py-5">
    <div class="row justify-content-center">
      <div class="col-lg-10 col-xl-8">

        <!-- 页面标题 -->
        <div class="text-center mb-6">
          <h1 class="display-4 fw-bold text-primary">分类</h1>
        </div>

        <!-- 使用 flex + gap 控制分类模块之间的间距 -->
        <div class="d-flex flex-column" style="gap: 3rem;">
          <section v-for="category in categories" :key="category.title">
            <!-- 一级分类标题 -->
            <h2 class="h3 fw-bold text-dark border-bottom pb-2 mb-4">
              {{ category.title }}
            </h2>

            <!-- 卡片网格 -->
            <div class="row g-4">
              <div v-for="item in category.items" :key="item.name" class="col-12 col-md-6 col-lg-4">
                <!-- 卡片：使用 flex-column + h-100 实现底部对齐 -->
                <div class="card h-100 shadow-sm border-0 d-flex flex-column">
                  <div class="card-body d-flex flex-column flex-grow-1">
                    <h3 class="h6 fw-bold text-dark mb-2">
                      {{ item.name }}
                    </h3>
                    <p class="text-muted fs-6 mb-3 flex-grow-1">
                      {{ item.desc }}
                    </p>
                    <!-- 按钮自动推到底部 -->
                    <a :href="item.url" target="_blank" rel="noopener" class="btn btn-outline-primary btn-sm mt-auto">
                      查看详情
                    </a>
                  </div>
                </div>
              </div>
            </div>
          </section>
        </div>

        <!-- 提示信息 -->
        <div class="text-center mt-5 pt-4 border-top">
          <p class="text-secondary fs-6">
            <i class="bi bi-info-circle me-2"></i>
            内容持续更新中，欢迎关注最新动态
          </p>
        </div>
      </div>
    </div>
  </div>
</template>

<script>
export default {
  name: 'CategoriesView',
  data() {
    return {
      categories: [
        {
          title: '笔记',
          items: [
            {
              name: '生物信息学基础',
              desc: '涵盖序列比对、数据库使用、基因组分析等核心知识点。',
              url: '/notes/bioinfo-basics'
            },
            {
              name: '结构生物学导论',
              desc: '蛋白质结构解析、PDB 数据分析、分子可视化工具实践。',
              url: '/notes/structural-bio'
            },
            {
              name: 'Python 科研脚本',
              desc: '常用数据处理、自动化脚本与生物信息流程开发技巧。',
              url: '/notes/python-for-research'
            }
          ]
        },
        {
          title: '项目',
          items: [
            {
              name: 'AlphaFold2 应用实践',
              desc: '使用 AlphaFold2 预测未知蛋白结构并进行功能分析。',
              url: '/projects/alphafold-practice'
            },
            {
              name: 'RNA-seq 差异分析流程',
              desc: '从原始数据到功能富集的全流程分析与可视化。',
              url: '/projects/rnaseq-pipeline'
            },
            {
              name: '分子对接模拟',
              desc: '基于 AutoDock Vina 的小分子-蛋白对接与结合模式分析。',
              url: '/projects/docking-study'
            }
          ]
        },
        {
          title: '课题',
          items: [
            {
              name: '膜蛋白结构与功能研究',
              desc: '聚焦 GPCR 与离子通道的结构动态与药物结合机制。',
              url: '/topics/membrane-protein-study'
            },
            {
              name: 'AI 在结构预测中的应用',
              desc: '探索深度学习模型在蛋白质折叠与设计中的前沿进展。',
              url: '/topics/ai-in-structural-bio'
            },
            {
              name: '进化与结构关系分析',
              desc: '从序列进化角度解析蛋白质结构保守性与变异规律。',
              url: '/topics/evolution-structure'
            }
          ]
        }
      ]
    }
  }
}
</script>

<style scoped></style>
