<!-- Resource.vue -->
<template>
  <div class="container py-4 px-3">
    <div class="row justify-content-center">
      <div class="col-lg-10 col-xl-8">
        <!-- 页面标题 -->
        <div class="text-center mb-4">
          <h1 class="display-6 fw-bold text-primary">资源</h1>
          <p class="text-muted mt-3 mb-0" style="font-size: 1.125rem; line-height: 1.5; color: #6c757d;">
            生物信息学与结构生物学领域常用工具
          </p>
        </div>

        <!-- 资源列表 -->
        <div class="d-flex flex-column" style="gap: 2rem;">
          <section v-for="category in resources" :key="category.title">
            <!-- 一级分类：浅灰下划线 -->
            <h2 class="h5 fw-semibold text-dark pb-2 mb-3" style="border-bottom: 1px solid #dee2e6;">
              {{ category.title }}
            </h2>

            <div class="ms-3">
              <div v-for="sub in category.children" :key="sub.title" class="mb-4">
                <!-- 二级分类 -->
                <h3 class="h6 fw-semibold text-dark mb-3">
                  {{ sub.title }}
                </h3>

                <!-- 资源条目 -->
                <ul class="list-unstyled mb-0">
                  <li v-for="item in sub.items" :key="item.name" class="mb-3">
                    <!-- 资源链接 -->
                    <a :href="item.url" target="_blank" rel="noopener noreferrer"
                      class="text-decoration-none text-primary fw-medium" style="font-size: 1rem;">
                      {{ item.name }}
                    </a>
                    <!-- 描述文字：标准博文大小 -->
                    <p class="text-muted mb-0 ms-3" style="font-size: 1rem; line-height: 1.6; color: #6c757d;">
                      {{ item.desc }}
                    </p>
                  </li>
                </ul>
              </div>
            </div>
          </section>
        </div>

        <!-- 提示信息：使用相同浅灰线 -->
        <div class="text-center mt-5 pt-4" style="border-top: 1px solid #dee2e6;">
          <p class="text-secondary mb-0" style="font-size: 0.875rem; color: #6c757d;">
            <i class="bi bi-info-circle me-1"></i>
            资源页面正在持续更新中，敬请期待更多优质内容
          </p>
        </div>
      </div>
    </div>
  </div>
</template>

<script>
import resources from '@/content/resources.json'

export default {
  name: 'ResourceView',
  data() {
    return { resources }
  }
}
</script>

<style scoped></style>
