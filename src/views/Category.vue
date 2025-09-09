<!-- Category.vue -->
<template>
  <div class="container py-4">
    <div class="row justify-content-center">
      <div class="col-lg-10 col-xl-8">
        <!-- 页面标题 -->
        <div class="text-center mb-4">
          <h1 class="display-6 fw-bold text-primary">{{ helper.pageTitle }}</h1>
        </div>

        <!-- 分类模块 -->
        <div class="d-flex flex-column" style="gap: 2rem;">
          <section v-for="(category, index) in categoryList" :key="index">
            <!-- 一级分类标题：浅灰下划线 -->
            <h2 class="h5 fw-semibold text-dark pb-2 mb-3" style="border-bottom: 1px solid #dee2e6;">
              {{ category.title }}
            </h2>

            <!-- 卡片网格 -->
            <div class="row g-3">
              <div v-for="(item, idx) in category.items" :key="idx" class="col-12 col-md-6 col-lg-4">
                <div class="card h-100 shadow-sm border-0 bg-white">
                  <div class="card-body p-4 d-flex flex-column">
                    <!-- 标题 -->
                    <h3 class="h6 fw-semibold text-dark mb-2">
                      {{ item.name }}
                    </h3>

                    <!-- 简介 -->
                    <p class="text-muted mb-3 flex-grow-1"
                      style="font-size: 0.875rem; line-height: 1.5; color: #6c757d;">
                      {{ item.desc }}
                    </p>

                    <!-- 分隔线：更浅的灰色 -->
                    <hr class="my-2" style="height: 1px; background-color: #dee2e6; opacity: 1; border: none;" />

                    <!-- 查看更多 -->
                    <div class="text-end">
                      <span class="text-primary fw-medium small d-inline-flex align-items-center gap-1 cursor-pointer"
                        :style="{ transition: 'transform 0.2s ease' }" @click="() => handleJump(item.url, item.type)"
                        @mouseenter="e => e.target.style.transform = 'translateX(2px)'"
                        @mouseleave="e => e.target.style.transform = 'translateX(0)'">
                        查看更多
                        <i class="bi bi-arrow-right"></i>
                      </span>
                    </div>
                  </div>
                </div>
              </div>
            </div>
          </section>
        </div>
      </div>
    </div>
  </div>
</template>

<script>
import categoryData from '@/content/categories.json';

export default {
  name: 'CategoryView',
  data() {
    return {
      categoryList: [],
      helper: {}
    };
  },
  created() {
    this.categoryList = categoryData.categoryList;
    this.helper = categoryData.categoryHelper;
  },
  methods: {
    handleJump(url, type) {
      if (!url) return;

      if (type === 'internal') {
        this.$router.push(url).catch(err => {
          if (err.name !== 'NavigationDuplicated' && !err.toString().includes('Navigation cancelled')) {
            console.error('Navigation error:', err);
          }
        });
      } else if (type === 'external') {
        window.open(url, '_blank', 'noopener,noreferrer');
      }
    }
  }
};
</script>

<style scoped></style>
