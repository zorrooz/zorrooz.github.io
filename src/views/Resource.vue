<!-- Resource.vue -->
<template>
  <div class="container py-5">
    <div class="row justify-content-center">
      <div class="col-lg-10 col-xl-8">
        <!-- 页面标题 -->
        <div class="text-center mb-6">
          <h1 class="display-4 fw-bold text-primary">资源</h1>
          <p class="lead fs-4 text-muted mt-3">
            生物信息学与结构生物学领域常用工具
          </p>
        </div>

        <!-- 资源列表 -->
        <section v-for="category in resources" :key="category.title" class="mb-6">
          <!-- 一级分类：黑色粗体 -->
          <h2 class="h3 fw-bold text-dark border-bottom pb-2 mb-4">
            {{ category.title }}
          </h2>

          <div class="ms-3">
            <div v-for="sub in category.children" :key="sub.title" class="mb-4">
              <!-- 二级分类：黑色粗体 -->
              <h3 class="h5 fw-bold text-dark mb-3">
                {{ sub.title }}
              </h3>

              <!-- 资源条目 -->
              <ul class="list-unstyled">
                <li v-for="item in sub.items" :key="item.name" class="mb-3">
                  <!-- 资源链接：蓝色，非粗体 -->
                  <a :href="item.url" target="_blank" rel="noopener" class="text-decoration-none fs-5 text-primary">
                    {{ item.name }}
                  </a>
                  <!-- 描述文字 -->
                  <p class="text-muted fs-6 ms-3 mb-0">
                    {{ item.desc }}
                  </p>
                </li>
              </ul>
            </div>
          </div>
        </section>

        <!-- 提示信息 -->
        <div class="text-center mt-5 pt-4 border-top">
          <p class="text-secondary fs-6">
            <i class="bi bi-info-circle me-2"></i>
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
